"""
Aaron's Archive v3.2 — DNA Analysis Engine
All core algorithms from scratch. No paid dependencies.
AI optional: Groq (gsk_...) or Gemini (AIza...) — both free tier.
BLAST uses NCBI free public API.
"""
from flask import Flask, request, jsonify, send_from_directory, Response
import math, re, json, urllib.request, urllib.parse, urllib.error, os, gzip, io, time
from collections import defaultdict, deque

app = Flask(__name__, static_folder=None)  # no static folder — prevents directory/source exposure

# ── Security ─────────────────────────────────────────────────────────────────
app.config['SECRET_KEY'] = os.environ.get('SECRET_KEY', os.urandom(32))
app.config['MAX_CONTENT_LENGTH'] = 2 * 1024 * 1024  # 2 MB max request
app.config['JSON_SORT_KEYS'] = False          # faster JSON serialization
app.config['JSONIFY_PRETTYPRINT_REGULAR'] = False  # compact JSON responses

# Compression is handled globally in the after_request hook below (see add_headers).

# Allowed domains for external fetches (SSRF protection)
_ALLOWED_FETCH_DOMAINS = {
    'blast.ncbi.nlm.nih.gov',
    'rest.uniprot.org',
    'alphafold.ebi.ac.uk',
    'search.rcsb.org',
    'files.rcsb.org',
    'api.groq.com',
    'generativelanguage.googleapis.com',
}

def _safe_fetch(url, timeout=30, **kwargs):
    """Only allow fetches to whitelisted domains. Accepts both URL strings and Request objects."""
    from urllib.parse import urlparse
    # urllib.request.Request objects have a .full_url attribute — extract it for domain check
    actual_url = url.full_url if isinstance(url, urllib.request.Request) else url
    domain = urlparse(actual_url).netloc.lower()
    if not any(domain == d or domain.endswith('.'+d) for d in _ALLOWED_FETCH_DOMAINS):
        raise ValueError(f'Fetch to {domain} not permitted')
    return urllib.request.urlopen(url, timeout=timeout, **kwargs)

# ── Simple in-memory rate limiter ────────────────────────────────────────────
# Protects expensive routes (BLAST, AI, external API calls) from abuse.
# Per-process only (fine for a single Railway/Render instance); resets on restart.
_rate_buckets = defaultdict(deque)
_RATE_LIMITS = {  # route_prefix: (max_requests, window_seconds)
    '/blast_submit':      (10, 60),
    '/structure_search':  (20, 60),
    '/phylo_tree':        (15, 60),
    '/ramachandran':      (20, 60),
    '/analyze':           (30, 60),
}

def _client_ip():
    # Vercel/Railway/Render all overwrite X-Forwarded-For at their edge and do
    # not forward client-supplied values, specifically to prevent IP spoofing —
    # so trusting the first value here is safe on those platforms. This
    # assumption does NOT hold if this app is ever deployed directly on a bare
    # VPS or behind a proxy that blindly forwards client headers; in that case
    # a client could inject a fake X-Forwarded-For to bypass rate limiting.
    fwd = request.headers.get('X-Forwarded-For', '')
    return fwd.split(',')[0].strip() if fwd else (request.remote_addr or 'unknown')

def _rate_limited(path):
    limit_conf = next((v for k, v in _RATE_LIMITS.items() if path.startswith(k)), None)
    if not limit_conf: return False
    max_req, window = limit_conf
    key = f'{_client_ip()}:{path}'
    now = time.time()
    bucket = _rate_buckets[key]
    while bucket and bucket[0] < now - window:
        bucket.popleft()
    if len(bucket) >= max_req:
        return True
    bucket.append(now)
    return False

@app.before_request
def _check_rate_limit():
    if _rate_limited(request.path):
        return jsonify({'error': 'Too many requests — please slow down and try again in a minute.'}), 429

def _body():
    """Safely extract the JSON request body as a dict, no matter what the
    client sends — malformed JSON, a JSON `null`, an empty body, or a JSON
    array/string instead of an object. Every route calls this instead of
    touching request.json directly, so a malformed request always gets a
    clean 400 from the route's own validation instead of an unhandled 500."""
    d = request.get_json(silent=True)
    return d if isinstance(d, dict) else {}

# ═══════════════════════════════════════════════
# CONSTANTS
# ═══════════════════════════════════════════════
CODON_TABLE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}
NN_PARAMS = {
    'AA':(-7.9,-22.2),'AT':(-7.2,-20.4),'AC':(-8.4,-22.4),'AG':(-7.8,-21.0),
    'TA':(-7.2,-21.3),'TT':(-7.9,-22.2),'TC':(-8.2,-22.2),'TG':(-8.5,-22.7),
    'CA':(-8.5,-22.7),'CT':(-7.8,-21.0),'CC':(-8.0,-19.9),'CG':(-10.6,-27.2),
    'GA':(-8.2,-22.2),'GT':(-8.4,-22.4),'GC':(-9.8,-24.4),'GG':(-8.0,-19.9),
}
RESTRICTION_ENZYMES = {
    'EcoRI':('GAATTC',1,'Cloning, common'),'BamHI':('GGATCC',1,'Cloning, common'),
    'HindIII':('AAGCTT',1,'Cloning, common'),'SalI':('GTCGAC',1,'Cloning'),
    'XhoI':('CTCGAG',1,'Cloning'),'NcoI':('CCATGG',1,'Expression vectors'),
    'NdeI':('CATATG',2,'Expression vectors'),'XbaI':('TCTAGA',1,'Cloning'),
    'SpeI':('ACTAGT',1,'Cloning'),'PstI':('CTGCAG',5,'Cloning'),
    'KpnI':('GGTACC',5,'Cloning'),'SacI':('GAGCTC',5,'Cloning'),
    'SmaI':('CCCGGG',3,'Blunt end'),'EcoRV':('GATATC',3,'Blunt end'),
    'ClaI':('ATCGAT',2,'Cloning'),'NotI':('GCGGCCGC',2,'8-cutter, rare'),
    'PacI':('TTAATTAA',5,'8-cutter, rare'),'AscI':('GGCGCGCC',2,'8-cutter, rare'),
    'FseI':('GGCCGGCC',6,'8-cutter, rare'),'AluI':('AGCT',2,'Frequent cutter'),
    'HaeIII':('GGCC',2,'Frequent cutter'),'MboI':('GATC',0,'Frequent cutter'),
    'TaqI':('TCGA',1,'Thermostable'),'MspI':('CCGG',1,'CpG analysis'),
    'HpaII':('CCGG',1,'CpG methylation'),'DraI':('TTTAAA',3,'AT-rich'),
    'PvuII':('CAGCTG',3,'Blunt end'),'ScaI':('AGTACT',3,'Blunt end'),
    'NheI':('GCTAGC',1,'Cloning'),'MluI':('ACGCGT',1,'Cloning'),
    'BglII':('AGATCT',1,'BamHI-compatible'),'SphI':('GCATGC',5,'Cloning'),
    'ApaI':('GGGCCC',5,'Cloning'),'NsiI':('ATGCAT',5,'Cloning'),
    'StuI':('AGGCCT',3,'Blunt end'),'SspI':('AATATT',3,'AT-rich'),
}
AA_MW = {'A':89.09,'R':174.20,'N':132.12,'D':133.10,'C':121.16,'E':147.13,'Q':146.15,
         'G':75.03,'H':155.16,'I':131.17,'L':131.17,'K':146.19,'M':149.21,'F':165.19,
         'P':115.13,'S':105.09,'T':119.12,'W':204.23,'Y':181.19,'V':117.15}
AA_PKA = {'D':(3.9,'neg'),'E':(4.1,'neg'),'H':(6.0,'pos'),'C':(8.3,'neg'),
          'Y':(10.1,'neg'),'K':(10.5,'pos'),'R':(12.5,'pos')}
AA_PKA_NT, AA_PKA_CT = 8.0, 3.1
AA_HYDRO = {'A':1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C':2.5,'E':-3.5,'Q':-3.5,'G':-0.4,
            'H':-3.2,'I':4.5,'L':3.8,'K':-3.9,'M':1.9,'F':2.8,'P':-1.6,'S':-0.8,
            'T':-0.7,'W':-0.9,'Y':-1.3,'V':4.2}
DIWV = {
    'AA':1.0,'AC':44.94,'AD':-7.49,'AE':1.0,'AF':1.0,'AG':1.0,'AH':-7.49,'AI':1.0,
    'AK':1.0,'AL':1.0,'AM':1.0,'AN':1.0,'AP':20.26,'AQ':1.0,'AR':1.0,'AS':1.0,
    'AT':1.0,'AV':1.0,'AW':1.0,'AY':1.0,'CA':1.0,'CC':1.0,'CD':20.26,'CE':1.0,
    'CF':1.0,'CG':1.0,'CH':33.6,'CI':1.0,'CK':1.0,'CL':20.26,'CM':33.6,'CN':1.0,
    'CP':20.26,'CQ':-6.54,'CR':1.0,'CS':1.0,'CT':33.6,'CV':-6.54,'CW':24.68,'CY':1.0,
    'DA':1.0,'DC':1.0,'DD':1.0,'DE':1.0,'DF':-6.54,'DG':1.0,'DH':1.0,'DI':1.0,
    'DK':-7.49,'DL':1.0,'DM':1.0,'DN':1.0,'DP':1.0,'DQ':1.0,'DR':-6.54,'DS':1.0,
    'DT':-14.03,'DV':1.0,'DW':1.0,'DY':1.0,'EA':1.0,'EC':44.94,'ED':20.26,'EE':33.6,
    'EF':1.0,'EG':1.0,'EH':-6.54,'EI':20.26,'EK':1.0,'EL':1.0,'EM':1.0,'EN':1.0,
    'EP':20.26,'EQ':20.26,'ER':1.0,'ES':20.26,'ET':1.0,'EV':1.0,'EW':-14.03,'EY':1.0,
    'FA':1.0,'FC':1.0,'FD':13.34,'FE':1.0,'FF':1.0,'FG':1.0,'FH':1.0,'FI':1.0,
    'FK':-14.03,'FL':1.0,'FM':1.0,'FN':1.0,'FP':20.26,'FQ':1.0,'FR':1.0,'FS':1.0,
    'FT':1.0,'FV':1.0,'FW':1.0,'FY':33.601,'GA':-7.49,'GC':1.0,'GD':1.0,'GE':-6.54,
    'GF':1.0,'GG':13.34,'GH':1.0,'GI':-7.49,'GK':-7.49,'GL':1.0,'GM':1.0,'GN':-7.49,
    'GP':1.0,'GQ':1.0,'GR':1.0,'GS':1.0,'GT':-7.49,'GV':1.0,'GW':13.34,'GY':-7.49,
    'HA':1.0,'HC':1.0,'HD':1.0,'HE':1.0,'HF':-9.37,'HG':-9.37,'HH':1.0,'HI':44.94,
    'HK':24.68,'HL':1.0,'HM':1.0,'HN':24.68,'HP':-1.88,'HQ':1.0,'HR':1.0,'HS':1.0,
    'HT':-6.54,'HV':1.0,'HW':-1.88,'HY':44.94,'IA':1.0,'IC':1.0,'ID':1.0,'IE':44.94,
    'IF':1.0,'IG':1.0,'IH':13.34,'II':1.0,'IK':-7.49,'IL':20.26,'IM':1.0,'IN':1.0,
    'IP':-1.88,'IQ':1.0,'IR':1.0,'IS':1.0,'IT':1.0,'IV':-7.49,'IW':1.0,'IY':1.0,
    'KA':1.0,'KC':1.0,'KD':1.0,'KE':1.0,'KF':1.0,'KG':-7.49,'KH':1.0,'KI':-7.49,
    'KK':1.0,'KL':-7.49,'KM':33.6,'KN':1.0,'KP':-6.54,'KQ':24.64,'KR':33.6,'KS':1.0,
    'KT':1.0,'KV':-7.49,'KW':1.0,'KY':1.0,'LA':1.0,'LC':1.0,'LD':1.0,'LE':1.0,
    'LF':1.0,'LG':1.0,'LH':1.0,'LI':1.0,'LK':-7.49,'LL':1.0,'LM':1.0,'LN':1.0,
    'LP':20.26,'LQ':33.6,'LR':20.26,'LS':1.0,'LT':1.0,'LV':1.0,'LW':24.68,'LY':1.0,
    'MA':13.34,'MC':1.0,'MD':1.0,'ME':1.0,'MF':1.0,'MG':1.0,'MH':58.28,'MI':1.0,
    'MK':1.0,'ML':1.0,'MM':-1.88,'MN':1.0,'MP':44.94,'MQ':-6.54,'MR':-6.54,'MS':44.94,
    'MT':-1.88,'MV':1.0,'MW':1.0,'MY':24.68,'NA':1.0,'NC':-1.88,'ND':1.0,'NE':1.0,
    'NF':-14.03,'NG':-14.03,'NH':1.0,'NI':44.94,'NK':24.68,'NL':1.0,'NM':1.0,'NN':1.0,
    'NP':-1.88,'NQ':-6.54,'NR':1.0,'NS':1.0,'NT':-7.49,'NV':1.0,'NW':-9.37,'NY':1.0,
    'PA':20.26,'PC':-6.54,'PD':-6.54,'PE':18.38,'PF':20.26,'PG':1.0,'PH':1.0,'PI':1.0,
    'PK':1.0,'PL':1.0,'PM':-6.54,'PN':1.0,'PP':20.26,'PQ':20.26,'PR':-6.54,'PS':20.26,
    'PT':1.0,'PV':20.26,'PW':-1.88,'PY':1.0,'QA':1.0,'QC':-6.54,'QD':20.26,'QE':20.26,
    'QF':-6.54,'QG':1.0,'QH':1.0,'QI':1.0,'QK':1.0,'QL':1.0,'QM':1.0,'QN':1.0,
    'QP':20.26,'QQ':20.26,'QR':1.0,'QS':44.94,'QT':1.0,'QV':-6.54,'QW':1.0,'QY':-6.54,
    'RA':1.0,'RC':1.0,'RD':1.0,'RE':1.0,'RF':1.0,'RG':-7.49,'RH':20.26,'RI':1.0,
    'RK':1.0,'RL':1.0,'RM':1.0,'RN':13.34,'RP':20.26,'RQ':20.26,'RR':58.28,'RS':44.94,
    'RT':1.0,'RV':1.0,'RW':58.28,'RY':-6.54,'SA':1.0,'SC':33.6,'SD':1.0,'SE':20.26,
    'SF':1.0,'SG':1.0,'SH':1.0,'SI':1.0,'SK':1.0,'SL':1.0,'SM':1.0,'SN':1.0,
    'SP':44.94,'SQ':20.26,'SR':20.26,'SS':20.26,'ST':1.0,'SV':1.0,'SW':1.0,'SY':1.0,
    'TA':1.0,'TC':1.0,'TD':1.0,'TE':20.26,'TF':13.34,'TG':-7.49,'TH':1.0,'TI':1.0,
    'TK':1.0,'TL':1.0,'TM':1.0,'TN':-14.03,'TP':1.0,'TQ':-6.54,'TR':1.0,'TS':1.0,
    'TT':1.0,'TV':1.0,'TW':-14.03,'TY':1.0,'VA':1.0,'VC':1.0,'VD':-14.03,'VE':1.0,
    'VF':1.0,'VG':-7.49,'VH':1.0,'VI':1.0,'VK':-1.88,'VL':1.0,'VM':1.0,'VN':1.0,
    'VP':20.26,'VQ':1.0,'VR':1.0,'VS':1.0,'VT':-7.49,'VV':1.0,'VW':1.0,'VY':-6.54,
    'WA':-14.03,'WC':1.0,'WD':1.0,'WE':1.0,'WF':1.0,'WG':-9.37,'WH':24.68,'WI':1.0,
    'WK':1.0,'WL':13.34,'WM':24.68,'WN':13.34,'WP':1.0,'WQ':1.0,'WR':1.0,'WS':1.0,
    'WT':-14.03,'WV':-7.49,'WW':1.0,'WY':1.0,'YA':1.0,'YC':1.0,'YD':1.0,'YE':-6.54,
    'YF':1.0,'YG':-7.49,'YH':13.34,'YI':1.0,'YK':1.0,'YL':1.0,'YM':44.94,'YN':1.0,
    'YP':13.34,'YQ':1.0,'YR':-15.91,'YS':1.0,'YT':-7.49,'YV':1.0,'YW':-9.37,'YY':13.34,
}
BLOSUM62 = {
    ('A','A'):4,('A','R'):-1,('A','N'):-2,('A','D'):-2,('A','C'):0,('A','Q'):-1,
    ('A','E'):-1,('A','G'):0,('A','H'):-2,('A','I'):-1,('A','L'):-1,('A','K'):-1,
    ('A','M'):-1,('A','F'):-2,('A','P'):-1,('A','S'):1,('A','T'):0,('A','W'):-3,
    ('A','Y'):-2,('A','V'):0,('R','R'):5,('R','N'):-1,('R','D'):-2,('R','C'):-3,
    ('R','Q'):1,('R','E'):0,('R','G'):-2,('R','H'):0,('R','I'):-3,('R','L'):-2,
    ('R','K'):2,('R','M'):-1,('R','F'):-3,('R','P'):-2,('R','S'):-1,('R','T'):-1,
    ('R','W'):-3,('R','Y'):-2,('R','V'):-3,('N','N'):6,('N','D'):1,('N','C'):-3,
    ('N','Q'):0,('N','E'):0,('N','G'):0,('N','H'):1,('N','I'):-3,('N','L'):-3,
    ('N','K'):0,('N','M'):-2,('N','F'):-3,('N','P'):-2,('N','S'):1,('N','T'):0,
    ('N','W'):-4,('N','Y'):-2,('N','V'):-3,('D','D'):6,('D','C'):-3,('D','Q'):0,
    ('D','E'):2,('D','G'):-1,('D','H'):-1,('D','I'):-3,('D','L'):-4,('D','K'):-1,
    ('D','M'):-3,('D','F'):-3,('D','P'):-1,('D','S'):0,('D','T'):-1,('D','W'):-4,
    ('D','Y'):-3,('D','V'):-3,('C','C'):9,('C','Q'):-3,('C','E'):-4,('C','G'):-3,
    ('C','H'):-3,('C','I'):-1,('C','L'):-1,('C','K'):-3,('C','M'):-1,('C','F'):-2,
    ('C','P'):-3,('C','S'):-1,('C','T'):-1,('C','W'):-2,('C','Y'):-2,('C','V'):-1,
    ('Q','Q'):5,('Q','E'):2,('Q','G'):-2,('Q','H'):0,('Q','I'):-3,('Q','L'):-2,
    ('Q','K'):1,('Q','M'):0,('Q','F'):-3,('Q','P'):-1,('Q','S'):0,('Q','T'):-1,
    ('Q','W'):-2,('Q','Y'):-1,('Q','V'):-2,('E','E'):5,('E','G'):-2,('E','H'):0,
    ('E','I'):-3,('E','L'):-3,('E','K'):1,('E','M'):-2,('E','F'):-3,('E','P'):-1,
    ('E','S'):0,('E','T'):-1,('E','W'):-3,('E','Y'):-2,('E','V'):-2,('G','G'):6,
    ('G','H'):-2,('G','I'):-4,('G','L'):-4,('G','K'):-2,('G','M'):-3,('G','F'):-3,
    ('G','P'):-2,('G','S'):0,('G','T'):-2,('G','W'):-2,('G','Y'):-3,('G','V'):-3,
    ('H','H'):8,('H','I'):-3,('H','L'):-3,('H','K'):-1,('H','M'):-2,('H','F'):-1,
    ('H','P'):-2,('H','S'):-1,('H','T'):-2,('H','W'):-2,('H','Y'):2,('H','V'):-3,
    ('I','I'):4,('I','L'):2,('I','K'):-1,('I','M'):1,('I','F'):0,('I','P'):-3,
    ('I','S'):-2,('I','T'):-1,('I','W'):-3,('I','Y'):-1,('I','V'):3,('L','L'):4,
    ('L','K'):-2,('L','M'):2,('L','F'):0,('L','P'):-3,('L','S'):-2,('L','T'):-1,
    ('L','W'):-2,('L','Y'):-1,('L','V'):1,('K','K'):5,('K','M'):-1,('K','F'):-3,
    ('K','P'):-1,('K','S'):0,('K','T'):-1,('K','W'):-3,('K','Y'):-2,('K','V'):-2,
    ('M','M'):5,('M','F'):0,('M','P'):-2,('M','S'):-1,('M','T'):-1,('M','W'):-1,
    ('M','Y'):-1,('M','V'):1,('F','F'):6,('F','P'):-4,('F','S'):-2,('F','T'):-2,
    ('F','W'):1,('F','Y'):3,('F','V'):-1,('P','P'):7,('P','S'):-1,('P','T'):-1,
    ('P','W'):-4,('P','Y'):-3,('P','V'):-2,('S','S'):4,('S','T'):1,('S','W'):-3,
    ('S','Y'):-2,('S','V'):-2,('T','T'):5,('T','W'):-2,('T','Y'):-2,('T','V'):0,
    ('W','W'):11,('W','Y'):2,('W','V'):-3,('Y','Y'):7,('Y','V'):-1,('V','V'):4,
}

# ═══════════════════════════════════════════════
# UTILITIES
# ═══════════════════════════════════════════════
def clean_seq(s): return re.sub(r'[^A-Za-z]','',s).upper()
# IUPAC nucleotide ambiguity codes (real sequencing data often contains these
# for uncertain base calls: R=A/G, Y=C/T, S=G/C, W=A/T, K=G/T, M=A/C,
# B=C/G/T, D=A/G/T, H=A/C/T, V=A/C/G, N=any)
_DNA_IUPAC = set('ATGCRYSWKMBDHVN')
_RNA_IUPAC = set('AUGCRYSWKMBDHVN')
_PROTEIN_AA = set('ACDEFGHIKLMNPQRSTVWYX*')
# Letters that ONLY ever mean an amino acid — they never appear in any DNA/RNA
# IUPAC code, so their presence rules out nucleic acid entirely.
_PROTEIN_EXCLUSIVE = _PROTEIN_AA - _DNA_IUPAC  # {E,F,I,L,P,Q,X,*}

def detect_type(s):
    if not s: return 'UNKNOWN'
    c = set(s)
    # No protein-exclusive letters present, and every character is a valid
    # nucleotide IUPAC code -> confidently DNA/RNA, even with ambiguity codes
    # or the occasional stray/typo character (previously a single unexpected
    # character like 'M' would wrongly flip a whole DNA sequence to PROTEIN).
    if not (c & _PROTEIN_EXCLUSIVE):
        if 'U' in c and 'T' not in c and c <= _RNA_IUPAC: return 'RNA'
        if c <= _DNA_IUPAC: return 'DNA'
    if c <= _PROTEIN_AA: return 'PROTEIN'
    return 'UNKNOWN'
def complement(s): return s.translate(str.maketrans('ATGCatgc','TACGtacg'))
def rev_complement(s): return complement(s)[::-1]
def dna_to_rna(s): return s.replace('T','U')
def translate(s):
    p=[]
    for i in range(0,len(s)-2,3):
        aa=CODON_TABLE.get(s[i:i+3],'?'); p.append(aa)
        if aa=='*': break
    return ''.join(p)

# ═══════════════════════════════════════════════
# DNA — BASIC STATS
# ═══════════════════════════════════════════════
def base_stats(s):
    n=len(s); cnt={b:s.count(b) for b in 'ATGCN'}
    gc=round((cnt['G']+cnt['C'])/n*100,2) if n else 0
    return {'length':n,'counts':cnt,'gc':gc,
            'at':round((cnt['A']+cnt['T'])/n*100,2) if n else 0,
            'purine':round((cnt['A']+cnt['G'])/n*100,2) if n else 0}

def molecular_weight_dna(s):
    w={'A':313.21,'T':304.19,'G':329.21,'C':289.18}
    return round(sum(w.get(b,0) for b in s)-61.96,1)

def melting_temp(s):
    n=len(s); a,t,g,c=s.count('A'),s.count('T'),s.count('G'),s.count('C')
    w=2*(a+t)+4*(g+c) if n<14 else 81.5+16.6*math.log10(0.05)+0.41*(g+c)/n*100-675/n
    dH=dS=0
    for i in range(n-1):
        if s[i:i+2] in NN_PARAMS: h,sv=NN_PARAMS[s[i:i+2]]; dH+=h; dS+=sv
    dH*=1000; dS+=-10.8; R=1.987; ct=250e-9
    nn=round((dH/(dS+R*math.log(ct/4)))-273.15,1) if dS else 0
    return {'wallace':round(w,1),'nearest_neighbor':nn}

def _adaptive_window(n, base=100, floor=20, target_points=18):
    # A fixed 100bp window leaves short sequences (mRNAs, small ORFs, etc.)
    # with only 0-2 sliding-window positions — too few for Chart.js to draw
    # a line at all. But shrinking the window too far the other direction
    # is just as broken: with only a handful of bases per window, G/C (or
    # A/T) counts saturate to the same value constantly by chance, so the
    # skew ratio pins at exactly ±1 for long stretches — a jagged square
    # wave instead of a real curve. Keep windows at least `floor` bases so
    # the ratio stays statistically meaningful, and only shrink below the
    # full 100bp when the sequence itself is short.
    if n <= base * 1.5:
        w = max(floor, n // target_points)
        return min(w, n) if n else floor
    return base

def gc_window(s,window=None):
    window = window or _adaptive_window(len(s))
    step=max(1,window//2)
    return [{'pos':i+window//2,'gc':round((s[i:i+window].count('G')+s[i:i+window].count('C'))/window*100,1)}
            for i in range(0,len(s)-window+1,step)]

def gc_skew_window(s,window=None):
    window = window or _adaptive_window(len(s))
    step=max(1,window//2); res=[]
    for i in range(0,len(s)-window+1,step):
        w=s[i:i+window]; g,c=w.count('G'),w.count('C')
        res.append({'pos':i+window//2,'skew':round((g-c)/(g+c),3) if (g+c) else 0})
    return res

def at_skew_window(s,window=None):
    window = window or _adaptive_window(len(s))
    step=max(1,window//2); res=[]
    for i in range(0,len(s)-window+1,step):
        w=s[i:i+window]; a,t=w.count('A'),(w.count('T')+w.count('U'))
        res.append({'pos':i+window//2,'skew':round((a-t)/(a+t),3) if (a+t) else 0})
    return res

def cumulative_gc_skew(s):
    cum=0; res=[]; step=max(1,len(s)//300)
    for i,b in enumerate(s):
        if b=='G': cum+=1
        elif b=='C': cum-=1
        if i%step==0: res.append({'pos':i,'skew':cum})
    return res

def dinucleotide_freq(s):
    all_di=[a+b for a in 'ATGC' for b in 'ATGC']; counts={}
    for i in range(len(s)-1):
        d=s[i:i+2]
        if all(b in 'ATGC' for b in d): counts[d]=counts.get(d,0)+1
    return {d:counts.get(d,0) for d in all_di}

def shannon_entropy(s):
    n=len(s); ent=-sum(p*math.log2(p) for b in set(s) if (p:=s.count(b)/n)>0)
    mx=math.log2(min(4,len(set(s)))) if len(set(s))>1 else 1
    return {'entropy':round(ent,3),'complexity_pct':round(ent/mx*100,1) if mx else 0}

# ═══════════════════════════════════════════════
# DNA — ORFs & TRANSLATION
# ═══════════════════════════════════════════════
def find_orfs(s,min_len=30):
    res=[]; stops={'TAA','TAG','TGA'}
    for ss,strand in [(s,'+'),(rev_complement(s),'-')]:
        for frame in range(3):
            i,start=frame,None
            while i+3<=len(ss):
                c=ss[i:i+3]
                if c=='ATG' and start is None: start=i
                elif c in stops and start is not None:
                    if i+3-start>=min_len:
                        orf=ss[start:i+3]; prot=translate(orf)
                        rs=start if strand=='+' else len(s)-(i+3)
                        re_=i+3 if strand=='+' else len(s)-start
                        res.append({'frame':frame+1,'strand':strand,'start':rs,'end':re_,
                                    'length':i+3-start,'protein_len':len(prot)-1,
                                    'protein':prot[:40]+('...' if len(prot)>40 else ''),
                                    'stop_codon':c})
                    start=None
                i+=3
    res.sort(key=lambda x:-x['length']); return res[:20]

def all_six_frames(s):
    rc=rev_complement(s); frames=[]
    for i in range(3):
        t=translate(s[i:]); frames.append({'frame':f'+{i+1}','translation':t[:80]+('…' if len(t)>80 else '')})
    for i in range(3):
        t=translate(rc[i:]); frames.append({'frame':f'-{i+1}','translation':t[:80]+('…' if len(t)>80 else '')})
    return frames

def codon_usage(s):
    cnts={}
    for i in range(0,len(s)-2,3):
        c=s[i:i+3]
        if c in CODON_TABLE: cnts[c]=cnts.get(c,0)+1
    aa_codons={}
    for codon,aa in CODON_TABLE.items(): aa_codons.setdefault(aa,[]).append(codon)
    rscu={}
    for codon,aa in CODON_TABLE.items():
        syns=aa_codons[aa]; total=sum(cnts.get(c,0) for c in syns); exp=total/len(syns) if syns else 0
        obs=cnts.get(codon,0)
        rscu[codon]={'count':obs,'aa':aa,'rscu':round(obs/exp,2) if exp else 0}
    return rscu

# ═══════════════════════════════════════════════
# DNA — STRUCTURAL FEATURES
# ═══════════════════════════════════════════════
def restriction_map(s):
    res=[]
    for enz,(pat,cut,note) in RESTRICTION_ENZYMES.items():
        p=pat.replace('N','[ATGC]').replace('W','[AT]').replace('R','[AG]').replace('Y','[CT]')
        sites=[m.start() for m in re.finditer(f'(?={p})',s)]
        if sites: res.append({'enzyme':enz,'pattern':pat,'sites':sites,'count':len(sites),'note':note,'cut_offset':cut})
    return sorted(res,key=lambda x:x['enzyme'])

def cpg_islands(s,window=200,step=50):
    islands,in_island,start=[],False,0
    for i in range(0,len(s)-window+1,step):
        w=s[i:i+window]; gc=(w.count('G')+w.count('C'))/window*100
        exp=(w.count('C')*w.count('G'))/window if window else 0
        oe=w.count('CG')/exp if exp else 0
        if gc>=50 and oe>=0.6:
            if not in_island: start,in_island=i,True
        else:
            if in_island: islands.append({'start':start,'end':i+window,'length':i+window-start}); in_island=False
    if in_island: islands.append({'start':start,'end':len(s),'length':len(s)-start})
    return islands

_COMP_CHAR={'A':'T','T':'A','G':'C','C':'G','a':'t','t':'a','g':'c','c':'g'}
def find_palindromes(s,min_len=4,max_len=8):
    """Find reverse-complement palindromes. Checks complementary character pairs
    with early exit on first mismatch, instead of constructing the full
    reverse-complement string for every candidate substring."""
    res=[]; n=len(s)
    for length in range(min_len,max_len+1,2):
        half=length//2
        for i in range(n-length+1):
            ok=True
            for k in range(half):
                if _COMP_CHAR.get(s[i+k]) != s[i+length-1-k]: ok=False; break
            if ok: res.append({'position':i,'sequence':s[i:i+length],'length':length})
    return res[:30]

def find_microsatellites(s,min_repeat=3):
    res=[]
    for ml in range(2,7):
        i=0
        while i<len(s)-ml:
            motif=s[i:i+ml]
            if len(set(motif))==1: i+=1; continue
            j,count=i+ml,1
            while j+ml<=len(s) and s[j:j+ml]==motif: count+=1; j+=ml
            if count>=min_repeat: res.append({'motif':motif,'start':i,'end':j,'repeats':count,'total_len':j-i,'type':f'({ml})'}); i=j
            else: i+=1
    res.sort(key=lambda x:-x['total_len']); return res[:20]

def find_repeats(s,min_len=8,max_len=30):
    """Find repeated (non-overlapping) subsequences of length min_len..max_len.
    Uses a single hash-based pass per length (O(n) per length) rather than a
    .find() scan per starting position (O(n) per position -> O(n^2) per length),
    which dominated total request time for longer sequences."""
    n=len(s); repeats,seen=[],set()
    for length in range(min_len,min(max_len,n//2)):
        first_seen={}
        for i in range(max(0,n-length)):
            sub=s[i:i+length]
            if sub in seen: continue
            j=first_seen.get(sub)
            if j is None: first_seen[sub]=i
            elif i-j>=length:
                repeats.append({'sequence':sub,'pos1':j,'pos2':i,'length':length}); seen.add(sub)
    repeats.sort(key=lambda x:-x['length']); return repeats[:15]

def find_promoter_elements(s):
    """Detect key regulatory DNA sequences."""
    els=[]
    patterns=[
        ('TATA box',r'TATA[AT]A[AT]','Core eukaryotic promoter. Found ~30bp before gene start. RNA polymerase anchors here.'),
        ('Kozak consensus',r'[AG]CC[AG]CCATGG','Strong translation start signal in eukaryotes. The ATG here is the start codon.'),
        ('Kozak-like',r'[AG].{2}ATG','Weaker translation start context around an ATG start codon.'),
        ('Shine-Dalgarno',r'AGGAGG|AAGGAG|GAGGAG','Prokaryotic ribosome binding site, ~10bp before start codon. Only in bacteria.'),
        ('CAAT box',r'CCAAT','Eukaryotic promoter element ~80bp upstream. Increases transcription frequency.'),
        ('GC box (Sp1)',r'GGGCGG','Binding site for Sp1 transcription factor. Common in housekeeping gene promoters.'),
        ('Splice donor',r'GT[GA]AGT','5\' end of an intron (non-coding insert). Tells the spliceosome where to cut.'),
        ('Splice acceptor',r'[CT]{6,}AG','3\' end of an intron. The AG dinucleotide is the actual splice point.'),
        ('AP-1 site',r'TGAC.CA','Binding site for AP-1 transcription factor. Responds to growth signals and stress.'),
        ('NF-kB site',r'GGGAC.TTCC','Binding site for NF-κB. Involved in immune response and inflammation.'),
    ]
    for name,pat,note in patterns:
        for m in re.finditer(pat,s):
            els.append({'type':name,'position':m.start(),'sequence':m.group()[:20],'note':note})
    els.sort(key=lambda x:x['position']); return els[:50]

# ═══════════════════════════════════════════════
# RNA SECONDARY STRUCTURE — Nussinov Algorithm
# ═══════════════════════════════════════════════
_RNA_PAIRS={('A','U'),('U','A'),('G','C'),('C','G'),('G','U'),('U','G')}
def can_pair(a,b):
    return (a,b) in _RNA_PAIRS

def rna_fold_nussinov(s,min_loop=3):
    """Predict RNA 2D structure using Nussinov dynamic programming."""
    rna=s.replace('T','U'); n=min(len(rna),120)
    rna=rna[:n]; dp=[[0]*n for _ in range(n)]
    for span in range(min_loop+1,n):
        for i in range(n-span):
            j=i+span
            best=max(dp[i+1][j],dp[i][j-1])
            if can_pair(rna[i],rna[j]):
                val=(dp[i+1][j-1] if i+1<=j-1 else 0)+1; best=max(best,val)
            for k in range(i+1,j): best=max(best,dp[i][k]+dp[k+1][j])
            dp[i][j]=best
    struct=['.']*n; stack=[(0,n-1)]
    while stack:
        i,j=stack.pop()
        if i>=j: continue
        if dp[i][j]==dp[i+1][j]: stack.append((i+1,j))
        elif dp[i][j]==dp[i][j-1]: stack.append((i,j-1))
        elif can_pair(rna[i],rna[j]) and dp[i][j]==(dp[i+1][j-1] if i+1<=j-1 else 0)+1:
            struct[i]='('; struct[j]=')'; stack.append((i+1,j-1))
        else:
            for k in range(i+1,j):
                if dp[i][j]==dp[i][k]+dp[k+1][j]: stack.append((i,k)); stack.append((k+1,j)); break
    db=''.join(struct); paired=db.count('(')
    return {'sequence':rna,'structure':db,'length':n,'base_pairs':paired,
            'mfe_approx':round(-paired*1.5,1),'unpaired':db.count('.')}

# ═══════════════════════════════════════════════
# PROTEIN ANALYSIS
# ═══════════════════════════════════════════════
def protein_mw(s):
    aa=[a for a in s if a in AA_MW and a!='*']
    return round(sum(AA_MW[a] for a in aa)-18.02*(len(aa)-1),1)

def isoelectric_point(s):
    aa=[a for a in s if a in AA_MW and a!='*']
    def charge(pH):
        c=1/(1+10**(pH-AA_PKA_NT))-1/(1+10**(AA_PKA_CT-pH))
        for a in aa:
            if a in AA_PKA:
                pka,sign=AA_PKA[a]
                c+=(1/(1+10**(pH-pka)) if sign=='pos' else -1/(1+10**(pka-pH)))
        return c
    lo,hi=0.0,14.0
    for _ in range(1000):
        mid=(lo+hi)/2
        if charge(mid)>0: lo=mid
        else: hi=mid
        if hi-lo<0.001: break
    return round((lo+hi)/2,2)

def gravy_score(s):
    vals=[AA_HYDRO[a] for a in s if a in AA_HYDRO]
    return round(sum(vals)/len(vals),3) if vals else 0

def aromaticity(s):
    return round(sum(s.count(a) for a in 'FYW')/len(s)*100,2) if s else 0

def extinction_coefficient(s):
    """Molar extinction coefficient at 280nm (Gill & von Hippel 1989 / Pace et al. 1995,
    the standard method used by ExPASy ProtParam). Assumes Tyr/Trp/Cys are the only
    UV-absorbing residues at 280nm."""
    nW,nY,nC = s.count('W'),s.count('Y'),s.count('C')
    reduced = nW*5500+nY*1490
    return {'reduced':reduced,'disulfide_bonds':reduced+(nC//2)*125,'nW':nW,'nY':nY,'nC':nC}

def aliphatic_index(s):
    """Ikai (1980) aliphatic index — relative volume occupied by aliphatic side chains
    (Ala,Val,Ile,Leu). Higher values indicate greater thermostability."""
    n=len(s)
    if not n: return 0
    return round(s.count('A')/n*100 + 2.9*s.count('V')/n*100 + 3.9*(s.count('I')+s.count('L'))/n*100,2)

def instability_index(s):
    sc=sum(DIWV.get(s[i]+s[i+1],1.0) for i in range(len(s)-1))
    return round(10/len(s)*sc,2) if s else 0

def aa_composition(s):
    total=len([a for a in s if a in AA_MW])
    return {a:{'count':s.count(a),'pct':round(s.count(a)/total*100,1) if total else 0}
            for a in 'ACDEFGHIKLMNPQRSTVWY'}

def secondary_structure_propensity(s):
    # Chou-Fasman (1978) conformational parameters. Verified against two
    # independent authoritative sources (Voet & Voet "Biochemistry" textbook
    # and Rockefeller University's prowl reference) after the previous table
    # was found to have significant transcription errors in several residues
    # (e.g. Met sheet propensity was 1.67 here vs the correct 1.05; Pro helix
    # propensity was 0.20 vs the correct 0.57).
    helix={'A':1.42,'L':1.21,'M':1.45,'E':1.51,'Q':1.11,'H':1.00,'K':1.16,'R':0.98,
           'V':1.06,'I':1.08,'C':0.70,'Y':0.69,'F':1.13,'W':1.08,'T':0.83,'S':0.77,
           'D':1.01,'N':0.67,'G':0.57,'P':0.57}
    sheet={'V':1.70,'I':1.60,'C':1.19,'Y':1.47,'F':1.38,'W':1.37,'L':1.30,'T':1.19,
           'M':1.05,'A':0.83,'R':0.93,'G':0.75,'D':0.54,'K':0.74,'S':0.75,'H':0.87,
           'N':0.89,'P':0.55,'E':0.37,'Q':1.10}
    h=round(sum(helix.get(a,1.0) for a in s)/len(s),3) if s else 0
    sh=round(sum(sheet.get(a,1.0) for a in s)/len(s),3) if s else 0
    return {'helix_score':h,'sheet_score':sh,'likely':'α-Helix dominant' if h>sh else 'β-Sheet dominant'}

def hydrophobicity_profile(s,window=9):
    """Kyte-Doolittle sliding window hydrophobicity, computed with a rolling
    sum (add the incoming residue, remove the outgoing one) instead of
    resumming the full window at every position."""
    n=len(s)
    if n<window: return []
    res=[]
    cur=sum(AA_HYDRO.get(a,0) for a in s[:window])
    res.append({'pos':window//2,'hydro':round(cur/window,3)})
    for i in range(1,n-window+1):
        cur+=AA_HYDRO.get(s[i+window-1],0)-AA_HYDRO.get(s[i-1],0)
        res.append({'pos':i+window//2,'hydro':round(cur/window,3)})
    return res

def charge_profile(s,window=7):
    """Sliding window net charge — positive=basic region, negative=acidic.
    Rolling sum, same approach as hydrophobicity_profile."""
    n=len(s)
    if n<window: return []
    charged={'K':1,'R':1,'H':0.5,'D':-1,'E':-1}
    res=[]
    cur=sum(charged.get(a,0) for a in s[:window])
    res.append({'pos':window//2,'charge':round(cur,1)})
    for i in range(1,n-window+1):
        cur+=charged.get(s[i+window-1],0)-charged.get(s[i-1],0)
        res.append({'pos':i+window//2,'charge':round(cur,1)})
    return res

def charge_vs_ph(s):
    """Return net charge of the whole protein at each pH 0.0-14.0 (step 0.2).
    Used for the Charge vs pH titration curve shown in the UI."""
    aa=[a for a in s if a in AA_MW and a!='*']
    def net_charge(pH):
        c=1/(1+10**(pH-AA_PKA_NT)) - 1/(1+10**(AA_PKA_CT-pH))
        for a in aa:
            if a in AA_PKA:
                pka,sign=AA_PKA[a]
                c+=(1/(1+10**(pH-pka)) if sign=='pos' else -1/(1+10**(pka-pH)))
        return round(c,2)
    return [{'ph':round(ph*0.2,1),'charge':net_charge(ph*0.2)} for ph in range(71)]


def transmembrane_helices(s,window=19,threshold=1.6):
    """Predict transmembrane helices by Kyte-Doolittle hydrophobicity."""
    prof=hydrophobicity_profile(s,window); helices,in_h,start=[],False,0
    for pt in prof:
        if pt['hydro']>=threshold:
            if not in_h: start,in_h=pt['pos'],True
        else:
            if in_h:
                length=pt['pos']-start
                if length>=15:
                    avg=round(sum(AA_HYDRO.get(s[i],0) for i in range(start,pt['pos']))/max(1,length),2)
                    helices.append({'start':start,'end':pt['pos'],'length':length,'avg_hydro':avg})
                in_h=False
    return helices

def signal_peptide(s):
    """Heuristic signal peptide prediction: N-region + hydrophobic H-region + C-region."""
    if len(s)<20: return {'detected':False,'score':0,'cleavage_site':None}
    n_reg=s[:5]; h_reg=s[5:min(25,len(s))]
    n_charged=sum(1 for a in n_reg if a in 'KRH')
    h_hydro=sum(AA_HYDRO.get(a,0) for a in h_reg)/len(h_reg) if h_reg else 0
    detected=h_hydro>1.6 and n_charged>=1; cleavage=None
    if detected:
        for i in range(15,min(35,len(s)-3)):
            if s[i] in 'AG' and s[i+1] in 'AGSTV': cleavage=i+1; break
    return {'detected':detected,'score':round(h_hydro+n_charged*0.5,2),
            'cleavage_site':cleavage,'h_region_hydro':round(h_hydro,2)}

# ═══════════════════════════════════════════════
# PRIMERS
# ═══════════════════════════════════════════════
def primer_design(s,l=20):
    fwd=s[:l]; rev=rev_complement(s[-l:])
    def sp(p):
        gc=(p.count('G')+p.count('C'))/len(p)*100; tm=melting_temp(p)['wallace']
        sc=0
        if 40<=gc<=60: sc+=3
        if 55<=tm<=65: sc+=3
        if p[-2:] in ('GC','CG','GG','CC'): sc+=2
        return sc,round(gc,1),round(tm,1)
    fs,fg,ft=sp(fwd); rs,rg,rt=sp(rev)
    return {'forward':{'sequence':fwd,'gc':fg,'tm':ft,'score':fs},
            'reverse':{'sequence':rev,'gc':rg,'tm':rt,'score':rs},
            'product_size':len(s)}

def check_hairpin(s,min_stem=4):
    n=len(s)
    for stem in range(min_stem,n//2+1):
        for i in range(n-2*stem-1):
            if s[i:i+stem]==rev_complement(s[n-stem:]): return True,stem
    return False,0

def check_self_dimer(s,min_m=4):
    rc=rev_complement(s)
    return any(rc[i:i+min_m] in s for i in range(len(s)-min_m+1))

def check_pair_dimer(f,r,min_m=4):
    rr=rev_complement(r)
    return any(f[i:i+min_m] in r or f[i:i+min_m] in rr for i in range(len(f)-min_m+1))

def advanced_primers(s,n_pairs=8,lengths=(18,20,22,24)):
    pairs=[]
    for fl in lengths:
        for rl in lengths:
            for fs in range(0,min(len(s)//3,80),4):
                fwd=s[fs:fs+fl]
                if len(fwd)<fl: continue
                for re_ in range(len(s),max(len(s)*2//3,len(s)-120),-4):
                    rs=re_-rl
                    if rs<fs+100: continue
                    rev=rev_complement(s[rs:re_])
                    if len(rev)<rl: continue
                    def sp(p):
                        gc=(p.count('G')+p.count('C'))/len(p)*100; tm=melting_temp(p)['nearest_neighbor']
                        sc=0
                        if 40<=gc<=60: sc+=3
                        if 55<=tm<=68: sc+=3
                        if p[-2:] in ('GC','CG','GG','CC'): sc+=2
                        if not any(p[i]*4 in p for i in range(len(p))): sc+=1
                        return sc,round(gc,1),round(tm,1)
                    fsc,fgc,ftm=sp(fwd); rsc,rgc,rtm=sp(rev)
                    fh,_=check_hairpin(fwd); rh,_=check_hairpin(rev)
                    fd=check_self_dimer(fwd); rd=check_self_dimer(rev); pd=check_pair_dimer(fwd,rev)
                    dtm=abs(ftm-rtm)
                    total=fsc+rsc-(2 if fh else 0)-(2 if rh else 0)-(2 if fd else 0)-(2 if rd else 0)-(3 if pd else 0)+(2 if dtm<=2 else 1 if dtm<=5 else 0)
                    pairs.append({'forward':{'sequence':fwd,'gc':fgc,'tm':ftm,'score':fsc,'hairpin':fh,'self_dimer':fd,'length':fl},
                                  'reverse':{'sequence':rev,'gc':rgc,'tm':rtm,'score':rsc,'hairpin':rh,'self_dimer':rd,'length':rl},
                                  'product_size':re_-fs,'delta_tm':round(dtm,1),'pair_dimer':pd,'total_score':total})
                    if len(pairs)>=200: break
                if len(pairs)>=200: break
            if len(pairs)>=200: break
    pairs.sort(key=lambda x:-x['total_score']); return pairs[:n_pairs]

# ═══════════════════════════════════════════════
# ALIGNMENT
# ═══════════════════════════════════════════════
def bscore(a,b): return BLOSUM62.get((a,b),BLOSUM62.get((b,a),-1))

def smith_waterman(s1,s2,st='DNA',gap=-8):
    s1,s2=s1[:300],s2[:300]; m,n=len(s1),len(s2)
    sf=(lambda a,b:2 if a==b else -1) if st=='DNA' else bscore
    H=[[0]*(n+1) for _ in range(m+1)]; best,bp=0,(0,0)
    for i in range(1,m+1):
        for j in range(1,n+1):
            H[i][j]=max(0,H[i-1][j-1]+sf(s1[i-1],s2[j-1]),H[i-1][j]+gap,H[i][j-1]+gap)
            if H[i][j]>best: best=H[i][j]; bp=(i,j)
    a1,a2,ma=[],[],[]; i,j=bp
    while i>0 and j>0 and H[i][j]>0:
        sc=sf(s1[i-1],s2[j-1])
        if H[i][j]==H[i-1][j-1]+sc: a1.append(s1[i-1]);a2.append(s2[j-1]);ma.append('|' if s1[i-1]==s2[j-1] else '.'); i-=1;j-=1
        elif H[i][j]==H[i-1][j]+gap: a1.append(s1[i-1]);a2.append('-');ma.append(' ');i-=1
        else: a1.append('-');a2.append(s2[j-1]);ma.append(' ');j-=1
    a1=''.join(reversed(a1));a2=''.join(reversed(a2));ms=''.join(reversed(ma))
    return {'score':best,'identity':round(ms.count('|')/max(len(ms),1)*100,1),
            'alignment_len':len(a1),'matches':ms.count('|'),'gaps':a1.count('-')+a2.count('-'),
            'seq1_aligned':a1[:120],'match_line':ms[:120],'seq2_aligned':a2[:120],'method':'Smith-Waterman (local)'}

def needleman_wunsch(s1,s2,st='DNA',gap=-2):
    s1,s2=s1[:200],s2[:200]; m,n=len(s1),len(s2)
    sf=(lambda a,b:2 if a==b else -1) if st=='DNA' else bscore
    dp=[[0]*(n+1) for _ in range(m+1)]
    for i in range(m+1): dp[i][0]=i*gap
    for j in range(n+1): dp[0][j]=j*gap
    for i in range(1,m+1):
        for j in range(1,n+1):
            dp[i][j]=max(dp[i-1][j-1]+sf(s1[i-1],s2[j-1]),dp[i-1][j]+gap,dp[i][j-1]+gap)
    a1,a2,ma=[],[],[]; i,j=m,n
    while i>0 or j>0:
        if i>0 and j>0 and dp[i][j]==dp[i-1][j-1]+sf(s1[i-1],s2[j-1]): a1.append(s1[i-1]);a2.append(s2[j-1]);ma.append('|' if s1[i-1]==s2[j-1] else '.'); i-=1;j-=1
        elif i>0 and dp[i][j]==dp[i-1][j]+gap: a1.append(s1[i-1]);a2.append('-');ma.append(' ');i-=1
        else: a1.append('-');a2.append(s2[j-1]);ma.append(' ');j-=1
    a1=''.join(reversed(a1));a2=''.join(reversed(a2));ms=''.join(reversed(ma))
    return {'score':dp[m][n],'identity':round(ms.count('|')/max(len(ms),1)*100,1),
            'alignment_len':len(a1),'matches':ms.count('|'),'gaps':a1.count('-')+a2.count('-'),
            'seq1_aligned':a1[:120],'match_line':ms[:120],'seq2_aligned':a2[:120],'method':'Needleman-Wunsch (global)'}

# ═══════════════════════════════════════════════
# BATCH FASTA
# ═══════════════════════════════════════════════
def parse_fasta(text):
    seqs,name,buf=[],None,[]
    for line in text.strip().split('\n'):
        line=line.strip()
        if not line: continue
        if line.startswith('>'):
            if name and buf: seqs.append({'name':name,'seq':''.join(buf)})
            name=line[1:].strip()[:60] or f'Seq_{len(seqs)+1}'; buf=[]
        else: buf.append(re.sub(r'[^A-Za-z]','',line))
    if name and buf: seqs.append({'name':name,'seq':''.join(buf)})
    return seqs

def batch_analyze(seqs):
    res=[]
    for e in seqs[:50]:
        raw=e['seq'].upper(); t=detect_type(raw)
        r={'name':e['name'],'type':t,'length':len(raw)}
        if t in ('DNA','RNA'):
            s=base_stats(raw)
            r.update({'gc':s['gc'],'at':s['at'],'mw':round(molecular_weight_dna(raw)/1000,2),
                      'tm_nn':melting_temp(raw)['nearest_neighbor'],'orfs':len(find_orfs(raw)),
                      'restriction_count':len(restriction_map(raw)),'entropy':shannon_entropy(raw)['complexity_pct']})
        elif t=='PROTEIN':
            r.update({'mw_kda':round(protein_mw(raw)/1000,3),'pi':isoelectric_point(raw),
                      'gravy':gravy_score(raw),'instability':instability_index(raw)})
        res.append(r)
    return res

# ═══════════════════════════════════════════════
# BLAST — NCBI free REST API
# ═══════════════════════════════════════════════
def _extract_json_object(text):
    """
    Extract first complete JSON object from text that may have surrounding HTML.
    Uses bracket-counter — avoids json.loads failing on trailing HTML.
    This is the KEY fix for NCBI BLAST responses.
    """
    start = text.find('{"BlastOutput2"')
    if start == -1: start = text.find('{')
    if start == -1: return None
    depth = 0; in_str = False; escape = False
    for i, c in enumerate(text[start:], start):
        if escape: escape = False; continue
        if c == '\\' and in_str: escape = True; continue
        if c == '"': in_str = not in_str; continue
        if in_str: continue
        if c == '{': depth += 1
        elif c == '}':
            depth -= 1
            if depth == 0:
                try: return json.loads(text[start:i+1])
                except: return None
    return None

def blast_submit(seq, program='blastn', database='nt'):
    url = 'https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi'
    seq_clean = re.sub(r'[^A-Za-z]','',seq)[:2000]
    params = {
        'CMD':'Put', 'PROGRAM':program, 'DATABASE':database,
        'QUERY':seq_clean,
        'HITLIST_SIZE':'20', 'EXPECT':'10',
        'TOOL':'aarons_archive', 'EMAIL':'noreply@aaronsarchive.app',
        'WORD_SIZE':'11' if program=='blastn' else '3'
    }
    data = urllib.parse.urlencode(params).encode()
    req = urllib.request.Request(url, data=data, headers={
        'Content-Type':'application/x-www-form-urlencoded',
        'User-Agent':'Mozilla/5.0 (compatible; AaronsArchive/3.2)'
    })
    try:
        with _safe_fetch(req, timeout=45) as r:
            raw = r.read().decode('utf-8','ignore')
    except Exception as e:
        raise ValueError(f'Cannot reach NCBI — check internet connection. ({e})')
    rid = re.search(r'RID\s*=\s*([A-Z0-9]+)', raw)
    rtoe = re.search(r'RTOE\s*=\s*(\d+)', raw)
    if not rid:
        # Try JSON form
        rid = re.search(r'"RID"\s*:\s*"([A-Z0-9]+)"', raw)
    if not rid:
        # Show first 400 chars for debugging
        preview = re.sub(r'<[^>]+>','',raw).strip()[:400]
        raise ValueError(f'NCBI returned no RID. Response: {preview}')
    # NCBI's RTOE is only an estimate, and it can occasionally come back
    # unusually large (seen: 451s) — most likely for unusual param
    # combinations, but there's no reason to sit and wait that long before
    # even the FIRST status check when most searches finish well under a
    # minute regardless. Cap it; the normal 8s polling loop takes over
    # after this first check either way, so nothing is lost by checking
    # sooner than NCBI's conservative suggestion.
    wait_raw = int(rtoe.group(1)) if rtoe else 30
    return {'rid':rid.group(1), 'wait':min(wait_raw, 60)}

def _txt(el, tag, default=''):
    """Safe XML element text extractor."""
    e = el.find(tag)
    return e.text.strip() if e is not None and e.text else default

_GENE_STOPWORDS = {'PROTEIN','ISOFORM','VARIANT','PRECURSOR','PARTIAL','UNCHARACTERIZED','PUTATIVE',
    'HYPOTHETICAL','LIKE','TYPE','FAMILY','MEMBER','SUBUNIT','COMPLEX','DOMAIN',
    'CONTAINING','RECEPTOR','FACTOR','HOMOLOG','HOMO','SAPIENS','MUS','MUSCULUS',
    'RATTUS','NORVEGICUS','PREDICTED','LOW','QUALITY','CHAIN','SUSCEPTIBILITY',
    'ANTIGEN','CELLULAR','TUMOR','BREAST','CANCER'}
_GENE_NAME_HINTS = [
    (re.compile(r'\bbreast cancer type 1\b', re.I), 'BRCA1'),
    (re.compile(r'\bbreast cancer type 2\b', re.I), 'BRCA2'),
    (re.compile(r'\bcellular tumor antigen p53\b', re.I), 'TP53'),
    (re.compile(r'\btumor protein p53\b', re.I), 'TP53'),
    (re.compile(r'\bepidermal growth factor receptor\b', re.I), 'EGFR'),
    (re.compile(r'\binsulin\b', re.I), 'INS'),
]
def _extract_gene_name(hit_def):
    """Best-effort gene symbol extraction from a BLAST hit definition line,
    for feeding the 'View 3D Structure' lookup. The original implementation
    only checked for a parenthesized symbol like "...(TP53)...", a format
    UniProt/SwissProt sometimes uses — but RefSeq (the default, recommended
    database) formats hit titles as "<description> [<organism>]" with no
    parentheses at all, e.g. "GTPase HRas [Homo sapiens]". That pattern
    never matched, so gene_name was silently empty for essentially every
    RefSeq hit and the 3D structure button never appeared. This tries,
    in order: SwissProt's GN= field, the legacy parenthesized form, a small
    set of well-known descriptive-name → symbol hints, then finally the
    word immediately preceding the organism bracket (handling "-like"
    suffixes), filtered against common non-gene description words.
    """
    m = re.search(r'\bGN=([A-Za-z0-9\-]{2,15})\b', hit_def)
    if m: return m.group(1).upper()
    m = re.search(r'\(([A-Z][A-Z0-9]{1,9})\)', hit_def)
    if m: return m.group(1).upper()
    for pat, sym in _GENE_NAME_HINTS:
        if pat.search(hit_def): return sym
    m = re.search(r'([A-Za-z][A-Za-z0-9]{1,9})(?:-like)?\s*\[', hit_def)
    if m:
        candidate = m.group(1).upper()
        if candidate not in _GENE_STOPWORDS and not candidate.isdigit():
            return candidate
    return ''

def _parse_blast_xml(xml_text, query_len_fallback=1):
    """Parse NCBI BLAST XML format into hits list.
    XML is plain text, never zipped, stable for 20+ years."""
    import xml.etree.ElementTree as ET
    # Strip DOCTYPE declaration which can confuse ElementTree
    xml_clean = re.sub(r'<!DOCTYPE[^>]+>', '', xml_text).strip()
    root = ET.fromstring(xml_clean)

    # Query length: try BlastOutput_query-len, else Iteration_query-len, else fallback
    query_len = int(_txt(root, 'BlastOutput_query-len') or query_len_fallback)
    hits = []

    for iteration in root.iter('Iteration'):
        qlen = int(_txt(iteration, 'Iteration_query-len') or query_len)
        for hit in iteration.findall('.//Hit'):
            hit_def   = _txt(hit, 'Hit_def', 'Unknown')
            accession = _txt(hit, 'Hit_accession', 'N/A')
            hsp = hit.find('.//Hsp')
            if hsp is None:
                continue
            align_len = max(int(_txt(hsp, 'Hsp_align-len', '1')), 1)
            identity  = int(_txt(hsp, 'Hsp_identity', '0'))
            q_from    = int(_txt(hsp, 'Hsp_query-from', '0'))
            q_to      = int(_txt(hsp, 'Hsp_query-to',   '0'))
            evalue    = _txt(hsp, 'Hsp_evalue', '1')
            try:
                ev = float(evalue)
                evalue_str = '0.0' if ev == 0 else (f'{ev:.2e}' if ev < 0.001 else f'{ev:.4f}')
            except:
                evalue_str = evalue
            gene_match = _extract_gene_name(hit_def)
            hits.append({
                'accession': accession,
                'title':     hit_def[:120],
                'gene_name': gene_match,
                'organism':  '',
                'taxid':     '',
                'score':     int(_txt(hsp, 'Hsp_score', '0')),
                'bit_score': round(float(_txt(hsp, 'Hsp_bit-score', '0')), 1),
                'evalue':    evalue_str,
                'identity':  round(identity / align_len * 100, 1),
                'align_len': align_len,
                'gaps':      int(_txt(hsp, 'Hsp_gaps', '0')),
                'coverage':  round(abs(q_to - q_from + 1) / qlen * 100, 1),
                'qseq':      _txt(hsp, 'Hsp_qseq')[:80],
                'midline':   _txt(hsp, 'Hsp_midline')[:80],
                'hseq':      _txt(hsp, 'Hsp_hseq')[:80],
            })
    return hits, max(query_len, 1)

def blast_poll(rid):
    """Poll NCBI BLAST using XML format.

    Previous versions gated the XML fetch behind a separate SearchInfo
    request, treating the search as "not ready" unless the literal
    substring "Status=READY" appeared in that response. That's brittle: if
    NCBI's status line ever differs by whitespace, casing, or format for a
    given request type, this would report WAITING forever -- indistinguishable
    from a genuinely slow search -- right up until our own timeout kills it.
    That's consistent with what's been observed: searches "timing out" at
    exactly our configured ceiling regardless of how fast NCBI actually
    finished. Removed that intermediate guess entirely. Every poll now goes
    straight for the real XML results, and "done" is decided the only way
    that can't be wrong: whether <BlastOutput> -- the actual results -- is
    present. If it isn't yet, we're still waiting, full stop, no format
    parsing involved.
    """
    url = 'https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi'
    _headers = {'User-Agent':'Mozilla/5.0 (compatible; AaronsArchive/3.2)','Accept':'text/html,application/xml'}

    try:
        req = urllib.request.Request(
            url + '?' + urllib.parse.urlencode({
                'CMD':'Get', 'FORMAT_TYPE':'XML', 'RID':rid,
                'FORMAT_OBJECT':'Alignment', 'HITLIST_SIZE':'20'
            }), headers=_headers)
        with _safe_fetch(req, timeout=60) as r:
            xml_raw = r.read().decode('utf-8','ignore')
    except urllib.error.HTTPError as e:
        return {'status':'FAILED','error':f'NCBI HTTP error {e.code}. Try again.'}
    except Exception as e:
        return {'status':'FAILED','error':f'Network error: {str(e)[:120]}'}

    if '<BlastOutput' in xml_raw:
        try:
            hits, query_len = _parse_blast_xml(xml_raw)
            return {'status':'DONE','hits':hits,'count':len(hits),'query_len':query_len}
        except Exception as e:
            return {'status':'FAILED','error':f'XML parse error: {type(e).__name__}: {str(e)[:150]}'}

    # Not ready yet -- check for a definitive failure/expiry signal using a
    # tolerant, case-insensitive, whitespace-flexible match (the exact
    # brittleness that caused this class of bug before). Anything else,
    # including a status line we don't recognize, safely defaults to still
    # running rather than guessing either way.
    if re.search(r'Status\s*=\s*FAILED', xml_raw, re.IGNORECASE):
        return {'status':'FAILED','error':'NCBI search failed on their servers'}
    if re.search(r'Status\s*=\s*UNKNOWN', xml_raw, re.IGNORECASE):
        return {'status':'UNKNOWN','error':'RID expired -- resubmit your sequence'}
    return {'status':'WAITING','debug':xml_raw[:200]}


def structure_search(gene_name, organism_id='9606'):
    """Search UniProt + AlphaFold + RCSB PDB for a gene's 3D structure.
    Returns structure metadata for the frontend viewer."""
    result = {'gene':gene_name,'uniprot_id':None,'uniprot_name':None,
              'alphafold_pdb_url':None,'alphafold_cif_url':None,
              'rcsb_ids':[],'function':None}
    try:
        # 1. UniProt lookup
        up_url = (f'https://rest.uniprot.org/uniprotkb/search'
                  f'?query=gene_exact:{urllib.parse.quote(gene_name)}'
                  f'+AND+organism_id:{organism_id}&format=json&size=1')
        req = urllib.request.Request(up_url,
            headers={'User-Agent':'AaronsArchive/3.2','Accept':'application/json'})
        with _safe_fetch(req, timeout=20) as r:
            up_data = json.loads(r.read().decode('utf-8','ignore'))
        if up_data.get('results'):
            entry = up_data['results'][0]
            uid = entry.get('primaryAccession','')
            result['uniprot_id'] = uid
            result['uniprot_name'] = entry.get('uniProtkbId','')
            # Extract function annotation
            for comment in entry.get('comments',[]):
                if comment.get('commentType') == 'FUNCTION':
                    texts = comment.get('texts',[])
                    if texts: result['function'] = texts[0].get('value','')[:300]; break
    except Exception as e:
        result['uniprot_error'] = str(e)[:80]

    if result['uniprot_id']:
        uid = result['uniprot_id']
        try:
            # 2. AlphaFold lookup
            af_url = f'https://alphafold.ebi.ac.uk/api/prediction/{uid}'
            req = urllib.request.Request(af_url,
                headers={'User-Agent':'AaronsArchive/3.2','Accept':'application/json'})
            with _safe_fetch(req, timeout=20) as r:
                af_data = json.loads(r.read().decode('utf-8','ignore'))
            if af_data and isinstance(af_data, list):
                result['alphafold_pdb_url'] = af_data[0].get('pdbUrl','')
                result['alphafold_cif_url'] = af_data[0].get('cifUrl','')
                result['alphafold_version'] = af_data[0].get('latestVersion','')
        except Exception as e:
            result['alphafold_error'] = str(e)[:80]

        try:
            # 3. RCSB PDB lookup for experimental structures
            rcsb_url = ('https://search.rcsb.org/rcsbsearch/v2/query'
                        '?json=' + urllib.parse.quote(json.dumps({
                'query':{'type':'terminal','service':'text',
                          'parameters':{'attribute':'rcsb_polymer_entity.rcsb_gene_name.value',
                                         'operator':'exact_match','value':gene_name}},
                'return_type':'entry',
                'request_options':{'pager':{'start':0,'rows':5},
                                     'sort':[{'sort_by':'score','direction':'desc'}]}
            })))
            req = urllib.request.Request(rcsb_url,
                headers={'User-Agent':'AaronsArchive/3.2','Accept':'application/json'})
            with _safe_fetch(req, timeout=15) as r:
                rcsb_data = json.loads(r.read().decode('utf-8','ignore'))
            result['rcsb_ids'] = [hit['identifier'] for hit in rcsb_data.get('result_set',[])[:5]]
        except Exception as e:
            result['rcsb_error'] = str(e)[:80]
    return result

# ═══════════════════════════════════════════════
# AI
# ═══════════════════════════════════════════════
def _build_prompt(stype, slen, r):
    p = [f"You are an expert molecular biologist and bioinformatician. Analyze this {stype} sequence ({slen} bp/aa) and provide a detailed, structured interpretation for a biotechnology student."]
    if stype in ('DNA', 'RNA'):
        p += [
            f"Key stats: GC content {r['stats']['gc']}%, AT content {r['stats']['at']}%.",
            f"Open Reading Frames detected: {len(r['orfs'])} (longest encodes {max((o['protein_len'] for o in r['orfs']), default=0)} aa).",
            f"Melting temperature (nearest-neighbor): {r['melting']['nearest_neighbor']}°C.",
            f"CpG islands: {len(r['cpg_islands'])}. Microsatellites/repeats: {len(r['microsatellites'])}.",
            f"Promoter elements found: {len(r.get('promoter_elements', []))}.",
            f"RNA secondary structure base pairs: {r.get('rna_fold', {}).get('base_pairs', 0)}.",
        ]
    else:
        p += [
            f"Molecular weight: {r['protein_mw']} Da. Isoelectric point (pI): {r['pi']}.",
            f"GRAVY score: {r['gravy']} ({'hydrophobic' if float(str(r['gravy']).replace(',','')) > 0 else 'hydrophilic'}).",
            f"Instability index: {r['instability']} ({'stable in vitro' if float(str(r['instability']).replace(',','')) < 40 else 'unstable in vitro'}).",
            f"Aliphatic index: {r.get('aliphatic_index','?')} (higher = more thermostable).",
            f"Extinction coefficient (280nm, reduced): {r.get('extinction',{}).get('reduced','?')} M^-1cm^-1.",
            f"Predicted transmembrane helices: {len(r.get('tm_helices', []))}.",
            f"Signal peptide: {'detected' if r.get('signal_peptide', {}).get('detected') else 'not detected'}.",
        ]
    p.append(
        "Provide a detailed analysis with the following sections:\n"
        "1. SEQUENCE IDENTITY — What type of gene/protein is this likely to be? What organism might it come from?\n"
        "2. FUNCTIONAL SIGNIFICANCE — What biological role does this sequence likely play?\n"
        "3. STRUCTURAL FEATURES — Comment on notable structural or compositional features.\n"
        "4. UNUSUAL FEATURES — Any anomalies, unusual GC content, regulatory signals, or interesting patterns?\n"
        "5. RESEARCH APPLICATIONS — How might a researcher use or study this sequence?\n"
        "Be specific, educational, and thorough. Use proper scientific terminology but explain key terms."
    )
    return '\n'.join(p)

def _groq(prompt, key):
    for model in ['llama-3.3-70b-versatile', 'llama-3.1-8b-instant', 'llama3-8b-8192']:
        try:
            req = urllib.request.Request('https://api.groq.com/openai/v1/chat/completions',
                data=json.dumps({'model': model, 'messages': [{'role': 'user', 'content': prompt}], 'max_tokens': 3000}).encode(),
                headers={'Content-Type': 'application/json', 'Authorization': f'Bearer {key}'})
            with _safe_fetch(req, timeout=30) as r:
                return {'text': json.loads(r.read())['choices'][0]['message']['content'], 'model': model}
        except urllib.error.HTTPError as e:
            b = e.read().decode('utf-8', 'ignore')
            if e.code == 401: return {'text': 'Invalid Groq key. Get one free at console.groq.com', 'model': ''}
            if e.code == 429: return {'text': 'Groq rate limit — wait a minute.', 'model': ''}
            if e.code != 404: return {'text': f'Groq error {e.code}: {b[:200]}', 'model': ''}
        except Exception as e:
            return {'text': f'Groq error: {e}', 'model': ''}
    return {'text': 'All Groq models unavailable.', 'model': ''}

def _gemini(prompt, key):
    # Current free-tier models (newest first)
    models = [
        'gemini-3.5-flash',
        'gemini-3.1-flash',
        'gemini-3.1-flash-lite',
        'gemini-3.1-pro',
        'gemini-2.5-flash',
        'gemini-2.5-pro',
        'gemini-2.0-flash',
    ]
    for model in models:
        try:
            url = f'https://generativelanguage.googleapis.com/v1beta/models/{model}:generateContent?key={key}'
            gen_config = {'maxOutputTokens': 4096, 'temperature': 0.7}
            # Gemini 2.5+/3.x models "think" before answering, silently eating the token
            # budget and truncating the visible answer. Disabling thinking (where supported)
            # ensures the full budget goes to the actual written response.
            if any(t in model for t in ('2.5', '3.')):
                gen_config['thinkingConfig'] = {'thinkingBudget': 0}
            body = json.dumps({'contents':[{'parts':[{'text':prompt}]}],
                                'generationConfig': gen_config}).encode()
            req = urllib.request.Request(url, data=body,
                headers={'Content-Type':'application/json','User-Agent':'AaronsArchive/3.2'})
            with _safe_fetch(req, timeout=30) as r:
                data = json.loads(r.read())
            candidates = data.get('candidates') or []
            if not candidates:
                fb = data.get('promptFeedback', {}).get('blockReason', 'no candidates returned')
                return {'text': f'Gemini returned no output (reason: {fb}). Try rephrasing or use a different model.', 'model': model}
            cand = candidates[0]
            finish = cand.get('finishReason', '')
            parts = cand.get('content', {}).get('parts', [])
            text = ''.join(p.get('text','') for p in parts if 'text' in p)
            if not text.strip():
                # Thinking-only response with no visible output — retry next model with thinking forced off
                if finish == 'MAX_TOKENS' and 'thinkingConfig' not in gen_config:
                    continue
                return {'text': f'Gemini ({model}) returned an empty response (finish reason: {finish}). Trying a different model may help.', 'model': model}
            note = ''
            if finish == 'MAX_TOKENS':
                note = '\n\n[Response was cut off at the token limit — the analysis above may be incomplete.]'
            return {'text': text + note, 'model': model}
        except urllib.error.HTTPError as e:
            b = e.read().decode('utf-8','ignore')
            if e.code == 400: return {'text': f'Google rejected this key: {b[:200]}', 'model': ''}
            if e.code == 403: return {'text': f'Gemini access denied. Details: {b[:200]}', 'model': ''}
            if e.code == 429: return {'text': 'Gemini rate limit hit — wait 1 minute and try again.', 'model': ''}
            if e.code == 404: continue  # model not found, try next
            return {'text': f'Gemini error {e.code}: {b[:200]}', 'model': ''}
        except Exception as e:
            return {'text': f'Gemini error: {str(e)[:150]}', 'model': ''}
    return {'text': 'No Gemini model available — try again later or use a Groq key (gsk_...) from console.groq.com', 'model': ''}

def _clean_api_key(raw):
    """Strip whitespace, invisible unicode characters, and common copy-paste
    artifacts (quotes, 'Bearer ' prefix, accidental labels) from a pasted API key."""
    key = raw.strip()
    # Strip zero-width/invisible unicode characters that .strip() doesn't catch
    for ch in ('\u200b', '\u200c', '\u200d', '\ufeff', '\u200e', '\u200f'):
        key = key.replace(ch, '')
    key = key.strip()
    # Strip surrounding quotes
    if len(key) >= 2 and key[0] == key[-1] and key[0] in ('"', "'"):
        key = key[1:-1].strip()
    # Strip accidental "Bearer " prefix
    if key.lower().startswith('bearer '):
        key = key[7:].strip()
    # Strip accidental label prefixes like "API Key:" or "Key:"
    key = re.sub(r'^(api[\s_-]?key|key)\s*[:=]\s*', '', key, flags=re.IGNORECASE)
    return key.strip()

def run_ai(prompt, key):
    key = _clean_api_key(key or '')
    if not key:
        return {'text': 'No API key provided.', 'model': 'none'}
    if len(key) < 8:
        return {'text': f'That doesn\'t look like a complete API key ({len(key)} characters). Double-check you copied the whole thing.', 'model': 'none'}
    if key.startswith('gsk_'):
        result = _groq(prompt, key)
        return result if isinstance(result, dict) else {'text': result, 'model': 'groq'}
    # Everything else is treated as a Gemini credential. We deliberately do NOT
    # hard-gate on a fixed prefix like 'AIza' here: Google has changed this format
    # before (rolling out 'AQ.'-prefixed Auth keys through 2026 as 'AIza' Standard
    # keys are phased out) and will likely do so again. Baking a specific prefix
    # into our own validation just means we reject perfectly valid keys the moment
    # the provider's format shifts. Instead we pass the key straight to Google's
    # API and let its own response be the source of truth on whether it's valid.
    result = _gemini(prompt, key)
    return result if isinstance(result, dict) else {'text': result, 'model': 'gemini'}

# ═══════════════════════════════════════════════
# PDF REPORT
# ═══════════════════════════════════════════════
def make_report(data):
    import datetime, html as _html
    now = datetime.datetime.now().strftime('%B %d, %Y  ·  %H:%M')
    t = data.get('type', '?')
    seq_preview = _html.escape((data.get('sequence_preview','') or data.get('sequence','') or '')[:180])

    css = """<style>
    @import url('https://fonts.googleapis.com/css2?family=Source+Serif+4:opsz,wght@8..60,500;8..60,600;8..60,700&family=IBM+Plex+Mono:wght@400;500;600;700&family=Inter:wght@400;500;600;700&display=swap');
    *{box-sizing:border-box;margin:0;padding:0}
    body{font-family:'Inter',Arial,sans-serif;font-size:11px;color:#22303c;background:#fff;-webkit-print-color-adjust:exact;print-color-adjust:exact}
    .page{max-width:960px;margin:0 auto;padding:0 0 40px}
    /* ── Letterhead ── */
    .hero{padding:30px 40px 20px;border-bottom:2.5px solid #0d3d52}
    .hero-top{display:flex;align-items:flex-end;justify-content:space-between}
    .hero-brand{display:flex;align-items:center;gap:12px}
    .hero-brand .mark{font-size:26px;line-height:1}
    .hero h1{font-family:'Source Serif 4',serif;font-size:23px;font-weight:600;color:#0d3d52;letter-spacing:0.002em}
    .hero .tagline{font-family:'Inter',sans-serif;font-size:9.5px;color:#7691a0;letter-spacing:0.08em;
      text-transform:uppercase;margin-top:2px;font-weight:500}
    .hero-badge{font-family:'IBM Plex Mono',monospace;font-size:9px;background:#eef6fa;
      border:1px solid #cfe7f0;color:#0077b6;padding:5px 13px;border-radius:3px;letter-spacing:0.05em;font-weight:600}
    .hero-meta{display:flex;gap:32px;margin-top:16px;flex-wrap:wrap;border-top:1px solid #e5edf1;padding-top:12px}
    .hero-meta div{font-family:'Inter',sans-serif;font-size:9px;color:#8398a5;text-transform:uppercase;letter-spacing:0.06em}
    .hero-meta b{color:#22303c;font-size:12px;display:block;margin-top:3px;font-weight:600;letter-spacing:0;text-transform:none;font-family:'IBM Plex Mono',monospace}
    .seq-preview{margin-top:14px;background:#f7fafb;border:1px solid #e5edf1;
      border-radius:5px;padding:9px 12px;font-family:'IBM Plex Mono',monospace;font-size:9px;color:#5c7385;
      letter-spacing:0.02em;word-break:break-all}
    /* ── Content ── */
    .content{padding:0 40px}
    h2{font-family:'Source Serif 4',serif;font-size:14.5px;color:#0d3d52;margin:28px 0 12px;
      font-weight:600;letter-spacing:0.001em;display:flex;align-items:center;gap:9px;
      border-bottom:1px solid #e5edf1;padding-bottom:8px}
    h2::before{content:'';width:4px;height:14px;background:#00b4d8;border-radius:1px;flex-shrink:0}
    h2 .cnt{font-family:'Inter',sans-serif;font-weight:600;color:#0077b6;
      background:#eef6fa;padding:1px 9px;border-radius:9px;font-size:10px;margin-left:auto}
    .summary{background:#f7fafb;border:1px solid #e5edf1;border-left:3px solid #0d3d52;
      border-radius:5px;padding:16px 19px;margin:8px 0 22px;line-height:1.9;font-size:10.8px;color:#3a4b57}
    .summary b{color:#0d3d52;font-weight:600}
    .summary ul{margin:7px 0 0 16px}
    .summary li{margin-top:5px}
    .summary-hd{font-family:'Inter',sans-serif;font-weight:700;color:#0d3d52;font-size:10.5px;
      text-transform:uppercase;letter-spacing:0.06em;display:block;margin-bottom:3px}
    .g{display:grid;grid-template-columns:repeat(4,1fr);gap:11px;margin-bottom:18px}
    .t{background:#fff;border:1px solid #e5edf1;border-radius:6px;padding:13px 10px;
      text-align:center;position:relative}
    .t::before{content:'';position:absolute;top:0;left:0;right:0;height:2.5px;background:#0d3d52;border-radius:6px 6px 0 0}
    .t .v{font-size:19px;font-weight:700;color:#0d3d52;font-family:'Source Serif 4',serif}
    .t .l{font-size:8.5px;color:#8398a5;text-transform:uppercase;letter-spacing:0.06em;margin-top:4px;font-weight:500}
    .compbar{display:flex;height:26px;border-radius:5px;overflow:hidden;margin:9px 0 5px;border:1px solid #e5edf1}
    .compbar div{display:flex;align-items:center;justify-content:center;font-family:'IBM Plex Mono',monospace;
      font-size:9px;color:#0a1620;font-weight:700}
    .complegend{display:flex;gap:18px;font-family:'Inter',sans-serif;font-size:9.5px;color:#5c7385;margin-bottom:18px}
    .complegend span{display:inline-flex;align-items:center;gap:5px}
    .complegend i{width:9px;height:9px;border-radius:2px;display:inline-block}
    table{width:100%;border-collapse:collapse;margin-bottom:18px;font-size:9.9px;
      border:1px solid #e5edf1;border-radius:6px;overflow:hidden}
    th{background:#0d3d52;color:#fff;padding:8px 10px;text-align:left;
      font-family:'Inter',sans-serif;font-size:8.8px;text-transform:uppercase;letter-spacing:0.05em;font-weight:600}
    td{padding:7px 10px;border-bottom:1px solid #eef3f6;font-family:'IBM Plex Mono',monospace;color:#3a4b57}
    tr:nth-child(even) td{background:#f9fbfc}
    tr:last-child td{border-bottom:none}
    code{background:#eef4f8;color:#0d3d52;padding:1.5px 5px;border-radius:3px;font-size:9.3px;font-family:'IBM Plex Mono',monospace}
    .ai{background:#fdfbf5;border:1px solid #ecdfb8;border-left:3px solid #c9971c;
      border-radius:5px;padding:16px 19px;margin-top:8px;line-height:1.9;font-size:10.6px;color:#5c4813}
    .ai-hd{font-family:'Inter',sans-serif;font-size:9.5px;color:#a17a10;text-transform:uppercase;
      letter-spacing:0.07em;margin-bottom:8px;font-weight:700}
    .footer{margin-top:38px;padding:16px 40px;border-top:1px solid #e5edf1;font-size:8.8px;color:#a3b2bb;
      display:flex;justify-content:space-between;font-family:'Inter',sans-serif}
    .methods{background:#fafbfc;border:1px solid #eef1f3;border-radius:5px;padding:14px 18px;margin-top:10px;
      font-size:9px;color:#7d8b95;line-height:1.8}
    .methods b{color:#5c6b76;font-weight:600}
    .methods-hd{font-family:'Inter',sans-serif;font-weight:700;color:#5c6b76;font-size:9.5px;
      text-transform:uppercase;letter-spacing:0.06em;display:block;margin-bottom:6px}
    .section{page-break-inside:avoid}
    .toc{display:flex;flex-wrap:wrap;gap:6px 18px;background:#f7fafb;border:1px solid #e5edf1;border-radius:6px;
      padding:12px 18px;margin:16px 0 4px;font-family:'IBM Plex Mono',monospace;font-size:9px;color:#5c7385}
    .toc b{color:#0d3d52;font-family:'Inter',sans-serif;text-transform:uppercase;letter-spacing:0.06em;
      font-size:8.5px;width:100%;margin-bottom:2px;font-weight:700}
    .toc span{color:#0077b6}
    .chart-img-box{background:#fff;border:1px solid #e5edf1;border-radius:8px;padding:12px;margin:6px 0 18px;
      page-break-inside:avoid}
    .chart-img-box img{width:100%;height:auto;display:block;border-radius:4px}
    .chart-img-cap{font-family:'Inter',sans-serif;font-size:8.5px;color:#8398a5;margin-top:8px;
      text-align:center;letter-spacing:0.02em}
    .chart-grid-2{display:grid;grid-template-columns:1fr 1fr;gap:14px}
    @page{margin:26px}
    </style>"""

    def stat_tile(v, l):
        return f'<div class="t"><div class="v">{v}</div><div class="l">{l}</div></div>'

    def table(headers, rows):
        h = ''.join(f'<th>{x}</th>' for x in headers)
        r = ''.join('<tr>' + ''.join(f'<td>{c}</td>' for c in row) + '</tr>' for row in rows)
        return f'<table><tr>{h}</tr>{r}</table>'

    charts = data.get('charts', {}) or {}
    def chart_img(chart_id, caption):
        src = charts.get(chart_id)
        if not src:
            return ''
        return f'<div class="chart-img-box"><img src="{src}" alt="{caption}"><div class="chart-img-cap">{caption}</div></div>'

    b = f'''<div class="page">
    <div class="hero">
      <div class="hero-top">
        <div class="hero-brand"><span class="mark">🧬</span><div><h1>Aaron's Archive</h1><div class="tagline">DNA · RNA · Protein Analysis Report</div></div></div>
        <div class="hero-badge">{_html.escape(str(t))} SEQUENCE</div>
      </div>
      <div class="hero-meta">
        <div>GENERATED<b>{now}</b></div>
        <div>LENGTH<b>{data.get("length",0):,} {"aa" if t=="PROTEIN" else "bp"}</b></div>
        <div>REPORT ID<b>AA-{abs(hash(str(data.get("sequence_preview","") or data.get("sequence",""))[:50])) % 100000:05d}</b></div>
      </div>
      {f'<div class="seq-preview">{seq_preview}{"…" if len(data.get("sequence","") or "")>180 else ""}</div>' if seq_preview else ''}
    </div>
    <div class="content">'''

    # ── Executive summary (rule-based, no fabricated claims) ──────────────────
    summary_pts = []
    if t in ('DNA','RNA'):
        s = data.get('stats',{}); m = data.get('melting',{})
        gc = s.get('gc',0)
        gc_desc = 'GC-rich' if gc>60 else 'AT-rich' if gc<40 else 'balanced'
        summary_pts.append(f'This {data.get("length",0):,} bp {t} sequence is <b>{gc_desc}</b> at <b>{gc}% GC</b>, with a nearest-neighbor melting temperature of <b>{m.get("nearest_neighbor","?")}°C</b>.')
        orfs = data.get('orfs',[])
        if orfs: summary_pts.append(f'<b>{len(orfs)} open reading frame(s)</b> were identified across all 6 frames, the longest encoding <b>{max((o.get("protein_len",0) for o in orfs), default=0)} aa</b>.')
        rest = data.get('restriction',[])
        if rest: summary_pts.append(f'<b>{sum(r.get("count",0) for r in rest)} restriction site(s)</b> found across <b>{len(rest)} enzyme(s)</b> screened.')
        cpg = data.get('cpg_islands',[])
        if cpg: summary_pts.append(f'<b>{len(cpg)} CpG island(s)</b> detected — often marking gene-regulatory regions.')
    elif t == 'PROTEIN':
        pi = data.get('pi',0); gravy = data.get('gravy',0)
        charge_desc = 'basic' if pi>7.5 else 'acidic' if pi<5.5 else 'near-neutral'
        hydro_desc = 'hydrophobic' if gravy>0.4 else 'hydrophilic' if gravy<-0.4 else 'moderately polar'
        summary_pts.append(f'This <b>{data.get("length",0):,} aa</b> protein has an isoelectric point of <b>pI {pi}</b> ({charge_desc}) and a GRAVY score of <b>{gravy}</b> ({hydro_desc}).')
        tm = data.get('tm_helices',[])
        if tm: summary_pts.append(f'<b>{len(tm)} transmembrane helix/helices</b> predicted — consistent with a membrane-associated protein.')
        if data.get('signal_peptide',{}).get('detected'): summary_pts.append('A <b>signal peptide</b> was detected, suggesting this protein is secreted or membrane-targeted.')
        inst = data.get('instability',0)
        summary_pts.append(f'Instability index of <b>{inst}</b> classifies this sequence as <b>{"unstable" if inst>40 else "stable"}</b> in vitro (Guruprasad scale).')

    if summary_pts:
        b += '<div class="summary"><span class="summary-hd">Summary</span><ul>' + ''.join(f'<li style="margin-top:4px">{p}</li>' for p in summary_pts) + '</ul></div>'

    # ── Table of contents (what's actually in this report) ────────────────
    toc_items = []
    if t in ('DNA','RNA'):
        toc_items.append('Core Statistics')
        if data.get('orfs'): toc_items.append(f"Open Reading Frames ({len(data['orfs'])})")
        if data.get('restriction'): toc_items.append(f"Restriction Sites ({len(data['restriction'])})")
        if data.get('promoter_elements'): toc_items.append(f"Regulatory Elements ({len(data['promoter_elements'])})")
        if data.get('cpg_islands'): toc_items.append(f"CpG Islands ({len(data['cpg_islands'])})")
        if data.get('microsatellites'): toc_items.append(f"Microsatellites ({len(data['microsatellites'])})")
        if data.get('palindromes'): toc_items.append(f"Palindromes ({len(data['palindromes'])})")
        if data.get('codon_usage'): toc_items.append('Codon Usage')
        if data.get('primers'): toc_items.append('Primer Design')
    elif t == 'PROTEIN':
        toc_items += ['Physicochemical Properties', 'Structure Overview']
        if data.get('tm_helices'): toc_items.append(f"Transmembrane Helices ({len(data['tm_helices'])})")
        if data.get('aa_composition'): toc_items.append('Amino Acid Composition')
    if data.get('ai_analysis'): toc_items.append('AI Interpretation')
    if toc_items:
        b += '<div class="toc"><b>Contents</b>' + ' <span>·</span> '.join(toc_items) + '</div>'

    if t in ('DNA','RNA'):
        s = data.get('stats',{}); m = data.get('melting',{})
        b += '<div class="section"><h2>Core Statistics</h2><div class="g">' + \
             stat_tile(f'{s.get("length",0):,}','Length (bp)') + \
             stat_tile(f'{s.get("gc",0)}%','GC Content') + \
             stat_tile(f'{m.get("nearest_neighbor",0)}°C','Tm (NN Method)') + \
             stat_tile(f'{data.get("mw",0)/1000:.2f} kDa','Mol. Weight') + '</div>' + \
             '<div class="g">' + \
             stat_tile(f'{m.get("wallace",0)}°C','Tm (Wallace)') + \
             stat_tile(f'{data.get("entropy",{}).get("complexity_pct","?")}%','Complexity') + \
             stat_tile(f'{len(data.get("microsatellites",[]))}','SSR Repeats') + \
             stat_tile(f'{len(data.get("palindromes",[]))}','Palindromes') + '</div>'

        counts = s.get('counts', {})
        total = sum(counts.get(k,0) for k in 'ATGC') or 1
        colors = {'A':'#4cf0b4','T':'#ff9494','G':'#ffc94d','C':'#52ecff'}
        bars = ''.join(f'<div style="width:{counts.get(k,0)/total*100:.2f}%;background:{colors[k]}">{k}: {counts.get(k,0)}</div>' for k in 'ATGC' if counts.get(k,0))
        legend = ''.join(f'<span><i style="background:{colors[k]}"></i>{k} — {counts.get(k,0)/total*100:.1f}%</span>' for k in 'ATGC')
        b += f'<div style="font-size:9.5px;color:#7691a0;font-family:\'IBM Plex Mono\',monospace;margin-bottom:4px">NUCLEOTIDE COMPOSITION</div><div class="compbar">{bars}</div><div class="complegend">{legend}</div></div>'
        b += chart_img('chartBase', 'Nucleotide composition breakdown')
        b += chart_img('chartGC', 'GC content — 100bp sliding window')

        orfs = data.get('orfs',[])
        if orfs:
            b += '<div class="section"><h2>Open Reading Frames<span class="cnt">' + str(len(orfs)) + '</span></h2>'
            b += table(['Frame','Strand','Start','End','Length (nt)','Protein (aa)'],
                       [[o["frame"],o["strand"],o["start"],o["end"],o["length"],o["protein_len"]] for o in orfs[:12]])
            b += '</div>'
        if charts.get('chartGCSkew') or charts.get('chartATSkew'):
            b += '<div class="section"><h2>Strand Skew Analysis</h2>'
            b += chart_img('chartGCSkew', 'GC skew — (G−C)/(G+C) per window')
            b += chart_img('chartATSkew', 'AT skew — (A−T)/(A+T) per window')
            b += chart_img('chartCumSkew', 'Cumulative GC skew — replication origin/terminus indicator')
            b += '</div>'
        rest = data.get('restriction',[])
        if rest:
            b += '<div class="section"><h2>Restriction Sites<span class="cnt">' + str(len(rest)) + '</span></h2>'
            b += table(['Enzyme','Recognition Site','Cuts','Note'],
                       [[f'<b>{r["enzyme"]}</b>', f'<code>{r["pattern"]}</code>', r["count"], r["note"]] for r in rest[:20]])
            b += '</div>'
        pe = data.get('promoter_elements',[])
        if pe:
            b += '<div class="section"><h2>Regulatory Elements<span class="cnt">' + str(len(pe)) + '</span></h2>'
            b += table(['Type','Position','Sequence','Note'],
                       [[e["type"], e["position"], f'<code>{e["sequence"]}</code>', e["note"]] for e in pe[:20]])
            b += '</div>'
        cpg = data.get('cpg_islands',[])
        if cpg:
            b += '<div class="section"><h2>CpG Islands<span class="cnt">' + str(len(cpg)) + '</span></h2>'
            b += table(['Start','End','Length','GC%','Obs/Exp'],
                       [[c.get("start"),c.get("end"),c.get("length"),c.get("gc"),c.get("obs_exp")] for c in cpg[:20]])
            b += '</div>'
        ssr = data.get('microsatellites',[])
        if ssr:
            b += '<div class="section"><h2>Microsatellites (SSRs)<span class="cnt">' + str(len(ssr)) + '</span></h2>'
            b += table(['Unit', 'Repeats', 'Position', 'Length'],
                       [[f'<code>{r.get("unit","")}</code>', r.get("count",r.get("repeats","")), r.get("position",r.get("pos","")), r.get("length","")] for r in ssr[:15]])
            b += '</div>'
        pal = data.get('palindromes',[])
        if pal:
            b += '<div class="section"><h2>Palindromic Sequences<span class="cnt">' + str(len(pal)) + '</span></h2>'
            b += table(['Sequence', 'Position', 'Length'],
                       [[f'<code>{p.get("sequence","")}</code>', p.get("position",p.get("pos","")), p.get("length","")] for p in pal[:15]])
            b += '</div>'
        codons = data.get('codon_usage',{})
        if codons:
            top_codons = sorted(((c,v) for c,v in codons.items() if v.get('aa') != '*'), key=lambda x:-x[1].get('count',0))[:12]
            if top_codons:
                b += '<div class="section"><h2>Codon Usage<span class="cnt">Top 12</span></h2>'
                b += table(['Codon','Amino Acid','Count','Frequency'],
                           [[f'<code>{c}</code>', v.get('aa','?'), v.get('count',0), f"{v.get('freq',0)}%" if 'freq' in v else '—'] for c,v in top_codons])
                b += chart_img('chartCodon', 'Top 20 most-used codons')
                b += chart_img('chartDi', 'Dinucleotide frequencies (CpG highlighted)')
                b += '</div>'
        primers = data.get('primers')
        if primers and isinstance(primers, dict) and primers.get('forward'):
            fwd, rev = primers.get('forward',{}), primers.get('reverse',{})
            b += '<div class="section"><h2>Primer Design</h2>'
            b += table(['Direction','Sequence','GC%','Tm (°C)'],
                       [['Forward', f'<code>{fwd.get("sequence","")}</code>', fwd.get('gc',''), fwd.get('tm','')],
                        ['Reverse', f'<code>{rev.get("sequence","")}</code>', rev.get('gc',''), rev.get('tm','')]])
            b += '</div>'
        if charts.get('chartMelt') or charts.get('chartEntropy'):
            b += '<div class="section"><h2>Thermodynamics &amp; Complexity</h2><div class="chart-grid-2">'
            b += chart_img('chartMelt', 'Simulated melting curve')
            b += chart_img('chartEntropy', 'Sequence complexity (Shannon entropy)')
            b += '</div></div>'

    elif t == 'PROTEIN':
        b += '<div class="section"><h2>Physicochemical Properties</h2><div class="g">' + \
             stat_tile(f'{data.get("protein_mw",0)/1000:.3f}','MW (kDa)') + \
             stat_tile(f'{data.get("pi",0)}','pI') + \
             stat_tile(f'{data.get("gravy",0)}','GRAVY') + \
             stat_tile(f'{data.get("instability",0)}','Instability') + '</div>' + \
             '<div class="g">' + \
             stat_tile(f'{data.get("aromaticity",0)}','Aromaticity') + \
             stat_tile(f'{data.get("aliphatic_index",0)}','Aliphatic Index') + \
             stat_tile(f'{data.get("extinction",{}).get("reduced","?") if isinstance(data.get("extinction"),dict) else data.get("extinction","?")}','Extinction (M⁻¹cm⁻¹)') + \
             stat_tile(f'{len(data.get("tm_helices",[]))}','TM Helices') + '</div></div>'
        sec = data.get('secondary',{})
        b += f'<div class="section"><h2>Structure Overview</h2><div class="summary" style="margin-top:0"><b>Likely fold:</b> {sec.get("likely","?")} &nbsp;·&nbsp; <b>TM helices:</b> {len(data.get("tm_helices",[]))} &nbsp;·&nbsp; <b>Signal peptide:</b> {"Yes" if data.get("signal_peptide",{}).get("detected") else "No"}</div>'
        b += chart_img('chartSecStr', 'Secondary structure propensity (α-helix vs β-sheet)')
        b += '</div>'
        tmh = data.get('tm_helices',[])
        if tmh:
            b += '<div class="section"><h2>Transmembrane Helices<span class="cnt">' + str(len(tmh)) + '</span></h2>'
            b += table(['#','Start','End','Length'], [[i+1,h.get("start"),h.get("end"),h.get("length")] for i,h in enumerate(tmh[:15])])
            b += '</div>'
        if charts.get('chartHydro'):
            b += '<div class="section"><h2>Hydrophobicity Profile</h2>'
            b += chart_img('chartHydro', 'Kyte-Doolittle hydrophobicity along the sequence (TM threshold = 1.6)')
            b += '</div>'
        if charts.get('chartCharge'):
            b += '<div class="section"><h2>Net Charge</h2>'
            b += chart_img('chartCharge', 'Net charge vs pH (or residue position)')
            b += '</div>'
        aac = data.get('aa_composition',{})
        if aac:
            top_aa = sorted(aac.items(), key=lambda x:-x[1].get('count',0))
            colors10 = ['#4cf0b4','#52ecff','#ffc94d','#dba6ff','#ff9494','#8fd9ee','#c9971c','#0077b6','#5c7385','#0d3d52']
            total_aa = sum(v.get('count',0) for _,v in top_aa) or 1
            bars = ''.join(f'<div style="width:{v.get("count",0)/total_aa*100:.2f}%;background:{colors10[i%10]}">{k if v.get("count",0)/total_aa>0.045 else ""}</div>' for i,(k,v) in enumerate(top_aa))
            b += '<div class="section"><h2>Amino Acid Composition</h2>' + f'<div class="compbar" style="height:22px">{bars}</div>'
            b += table(['Residue','Count','%'], [[k, v.get('count',0), f"{v.get('pct',round(v.get('count',0)/total_aa*100,1))}%"] for k,v in top_aa[:20]])
            b += chart_img('chartAA', 'Amino acid composition (%)')
            b += '</div>'

    ai = data.get('ai_analysis')
    if ai:
        ai_safe = _html.escape(str(ai)).replace('\n','<br>')
        b += f'<div class="section"><h2>AI Interpretation</h2><div class="ai"><div class="ai-hd">AI-Generated — verify independently</div>{ai_safe}</div></div>'

    # ── Methods note — grounds every figure in a named, checkable method ──
    if t in ('DNA','RNA'):
        b += '''<div class="methods"><span class="methods-hd">Methods</span>
        <b>Tm (NN):</b> nearest-neighbor thermodynamic method. &nbsp;<b>Tm (Wallace):</b> Wallace rule, 4(G+C)+2(A+T), for short primers.
        &nbsp;<b>GC skew:</b> (G−C)/(G+C) per window. &nbsp;<b>Complexity:</b> Shannon entropy of base composition, normalized to 0–100%.
        &nbsp;<b>ORFs:</b> ATG→stop, all 6 reading frames, minimum length filtered. &nbsp;<b>CpG islands:</b> Gardiner-Garden &amp; Frommer criteria (GC% &gt; 50, Obs/Exp CpG &gt; 0.6, ≥200bp).</div>'''
    elif t == 'PROTEIN':
        b += '''<div class="methods"><span class="methods-hd">Methods</span>
        <b>pI:</b> Henderson-Hasselbalch iterative estimation over standard pKa values. &nbsp;<b>GRAVY:</b> Kyte-Doolittle hydropathy, mean per residue.
        &nbsp;<b>Instability index:</b> Guruprasad et al. 1990 dipeptide method (&gt;40 = predicted unstable in vitro).
        &nbsp;<b>Aliphatic index:</b> Ikai 1980, relative volume of aliphatic side chains. &nbsp;<b>TM helices:</b> hydrophobicity window scan, threshold 1.6.</div>'''

    b += '</div><div class="footer"><span>Aaron\'s Archive — DNA Analysis Engine</span><span>' + now + '</span></div></div>'
    return f"<!DOCTYPE html><html><head><meta charset='UTF-8'><title>Aaron's Archive Report</title>{css}</head><body>{b}</body></html>"

# ═══════════════════════════════════════════════
# PHYLOGENETIC TREE — pure Python UPGMA
# ═══════════════════════════════════════════════
def _seq_distance(a, b):
    """Simple p-distance between two aligned sequences (proportion of differing sites)."""
    a = re.sub(r'[^A-Za-z]','',a).upper()
    b = re.sub(r'[^A-Za-z]','',b).upper()
    if not a or not b: return 1.0
    # Align by padding shorter sequence
    length = max(len(a),len(b))
    a = a.ljust(length,'-'); b = b.ljust(length,'-')
    diffs = sum(1 for x,y in zip(a,b) if x!=y and x!='-' and y!='-')
    comp  = sum(1 for x,y in zip(a,b) if x!='-' and y!='-')
    return diffs/comp if comp>0 else 1.0

def build_phylo_tree(sequences):
    """Build UPGMA phylogenetic tree from a dict of {name: sequence}.
    Returns Newick string + distance matrix for frontend rendering."""
    names = list(sequences.keys())
    seqs  = [sequences[n] for n in names]
    n     = len(names)
    if n < 2: raise ValueError('Need at least 2 sequences')
    if n > 50: raise ValueError('Max 50 sequences for tree building')

    # Build distance matrix
    dm = [[0.0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1,n):
            d = _seq_distance(seqs[i], seqs[j])
            dm[i][j] = dm[j][i] = d

    # UPGMA clustering
    clusters  = {i: [i] for i in range(n)}
    labels    = {i: names[i] for i in range(n)}
    newick_s  = {i: names[i].replace('(','_').replace(')','_').replace(',','_').replace(':','_') for i in range(n)}
    heights   = {i: 0.0 for i in range(n)}
    dist      = {}
    for i in range(n):
        for j in range(i+1,n):
            dist[(i,j)] = dm[i][j]
    active    = list(range(n))
    node_id   = n
    steps     = []  # for frontend rendering

    while len(active) > 1:
        min_d = float('inf'); pair = None
        for i in range(len(active)):
            for j in range(i+1,len(active)):
                a2,b2 = active[i],active[j]
                k = (min(a2,b2),max(a2,b2))
                if dist.get(k, float('inf')) < min_d:
                    min_d = dist[k]; pair=(a2,b2)
        a,b = pair
        h = min_d/2
        branch_a = round(h - heights[a], 5)
        branch_b = round(h - heights[b], 5)
        newick_s[node_id] = f'({newick_s[a]}:{branch_a},{newick_s[b]}:{branch_b})'
        labels[node_id]   = f'Node{node_id}'
        heights[node_id]  = h
        steps.append({'a':labels[a],'b':labels[b],'dist':round(min_d,4),'height':round(h,4)})
        
        sa = len(clusters[a]); sb = len(clusters[b])
        clusters[node_id] = clusters[a]+clusters[b]
        active.remove(a); active.remove(b)
        for c in active:
            da = dist.get((min(a,c),max(a,c)),0)
            db = dist.get((min(b,c),max(b,c)),0)
            dist[(min(node_id,c),max(node_id,c))] = (sa*da+sb*db)/(sa+sb)
        active.append(node_id)
        node_id += 1

    root = active[0]
    newick = newick_s[root]+';'
    return {
        'newick': newick,
        'names': names,
        'distance_matrix': dm,
        'steps': steps,
        'n': n
    }

# ═══════════════════════════════════════════════
# RAMACHANDRAN — phi/psi from PDB backbone
# ═══════════════════════════════════════════════
def _calc_dihedral(p1,p2,p3,p4):
    import math
    def sub(a,b): return [a[i]-b[i] for i in range(3)]
    def cross(a,b): return [a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]
    def dot(a,b):   return sum(a[i]*b[i] for i in range(3))
    def norm(a):
        m=math.sqrt(sum(x*x for x in a))
        return [x/m for x in a] if m>0 else a
    b1=sub(p2,p1); b2=sub(p3,p2); b3=sub(p4,p3)
    n1=cross(b1,b2); n2=cross(b2,b3)
    m1=cross(n1,norm(b2))
    return math.degrees(math.atan2(dot(m1,n2), dot(n1,n2)))

def ramachandran_from_pdb(pdb_url):
    """Fetch a PDB file and calculate phi/psi backbone dihedral angles.
    Returns list of {res, phi, psi, aa} for Ramachandran plot."""
    req = urllib.request.Request(pdb_url,
        headers={'User-Agent':'AaronsArchive/3.2','Accept':'text/plain'})
    with _safe_fetch(req, timeout=30) as r:
        pdb_text = r.read().decode('utf-8','ignore')

    # Parse ATOM records for backbone N, CA, C atoms (chain A, first model only)
    residues = {}  # resnum → {N,CA,C: [x,y,z]}
    in_model1 = True
    for line in pdb_text.splitlines():
        if line.startswith('ENDMDL'): break   # only first model
        if not line.startswith('ATOM'): continue
        try:
            atom_name = line[12:16].strip()
            chain     = line[21]
            if chain not in ('A',' '): continue
            resnum    = int(line[22:26].strip())
            res_name  = line[17:20].strip()
            x,y,z     = float(line[30:38]),float(line[38:46]),float(line[46:54])
            if atom_name in ('N','CA','C'):
                if resnum not in residues:
                    residues[resnum] = {'name':res_name}
                residues[resnum][atom_name] = [x,y,z]
        except: continue

    # Calculate phi/psi for each residue
    aa_1letter = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E',
                  'GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F',
                  'PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
    resnums = sorted(residues.keys())
    points  = []
    for idx, rn in enumerate(resnums):
        res = residues[rn]
        if not all(k in res for k in ('N','CA','C')): continue
        aa3  = res.get('name','UNK')
        aa1  = aa_1letter.get(aa3, 'X')
        phi = psi = None
        # Phi = dihedral C(i-1)-N(i)-CA(i)-C(i)
        if idx > 0:
            prev = residues[resnums[idx-1]]
            if 'C' in prev:
                try: phi = round(-_calc_dihedral(prev['C'],res['N'],res['CA'],res['C']),1)
                except: pass
        # Psi = dihedral N(i)-CA(i)-C(i)-N(i+1)
        if idx < len(resnums)-1:
            nxt = residues[resnums[idx+1]]
            if 'N' in nxt:
                try: psi = round(-_calc_dihedral(res['N'],res['CA'],res['C'],nxt['N']),1)
                except: pass
        if phi is not None and psi is not None:
            points.append({'res':rn,'aa':aa1,'aa3':aa3,'phi':phi,'psi':psi})

    if not points: raise ValueError('No backbone atoms found in PDB')
    return points

# ═══════════════════════════════════════════════
# ROUTES
# ═══════════════════════════════════════════════
@app.route('/')
def index():
    resp = send_from_directory('.','index.html')
    resp.headers['Cache-Control'] = 'public, max-age=300'  # cache 5 min
    return resp

@app.route('/analyze',methods=['POST'])
def analyze():
    d=_body(); raw=d.get('sequence','').strip()
    if not raw: return jsonify({'error':'No sequence provided'}),400
    seq=clean_seq(raw)
    if len(seq)>150000: return jsonify({'error':'Sequence too long (max 150,000 characters). Use Batch FASTA for multiple large sequences.'}),400
    if not seq: return jsonify({'error':'No valid sequence characters found'}),400
    st=detect_type(seq)
    res={'type':st,'length':len(seq),'sequence_preview':seq[:100]}
    if st in ('DNA','RNA'):
        if st=='RNA': seq=seq.replace('U','T')
        res.update({'stats':base_stats(seq),'mw':molecular_weight_dna(seq),'melting':melting_temp(seq),
            'complement':complement(seq)[:100],'rev_complement':rev_complement(seq)[:100],
            'six_frames':all_six_frames(seq),'orfs':find_orfs(seq),'restriction':restriction_map(seq),
            'cpg_islands':cpg_islands(seq),'palindromes':find_palindromes(seq),
            'microsatellites':find_microsatellites(seq),'repeats':find_repeats(seq),
            'entropy':shannon_entropy(seq),'dinucleotide':dinucleotide_freq(seq),
            'codon_usage':codon_usage(seq),'gc_window':gc_window(seq),
            'gc_skew':gc_skew_window(seq),'at_skew':at_skew_window(seq),
            'cumulative_skew':cumulative_gc_skew(seq),'promoter_elements':find_promoter_elements(seq),
            'rna_fold':rna_fold_nussinov(seq),'primers':primer_design(seq)})
        if st=='RNA': res['rna_sequence']=dna_to_rna(seq)[:100]
    elif st=='PROTEIN':
        res.update({'protein_mw':protein_mw(seq),'pi':isoelectric_point(seq),'gravy':gravy_score(seq),
            'aromaticity':aromaticity(seq),'instability':instability_index(seq),
            'extinction':extinction_coefficient(seq),'aliphatic_index':aliphatic_index(seq),
            'aa_composition':aa_composition(seq),'secondary':secondary_structure_propensity(seq),
            'hydro_profile':hydrophobicity_profile(seq),'charge_profile':charge_profile(seq),'charge_vs_ph':charge_vs_ph(seq),
            'tm_helices':transmembrane_helices(seq),'signal_peptide':signal_peptide(seq),
            'charge_at_7':round(sum(1/(1+10**(7-AA_PKA[a][0])) if AA_PKA[a][1]=='pos' else -1/(1+10**(AA_PKA[a][0]-7)) for a in seq if a in AA_PKA),2)})
    ai_key=d.get('api_key','').strip()
    ai_result = run_ai(_build_prompt(st,len(seq),res),ai_key) if ai_key else None
    res['ai_analysis'] = ai_result.get('text','') if isinstance(ai_result,dict) else ai_result
    res['ai_model']    = ai_result.get('model','') if isinstance(ai_result,dict) else ''
    return jsonify(res)

@app.route('/batch',methods=['POST'])
def batch():
    d=_body(); raw=d.get('fasta','').strip()
    if not raw: return jsonify({'error':'No FASTA data'}),400
    seqs=parse_fasta(raw)
    if not seqs: return jsonify({'error':'No valid sequences. Format: >name\\nSEQUENCE'}),400
    return jsonify({'count':len(seqs),'results':batch_analyze(seqs)})

@app.route('/align',methods=['POST'])
def align():
    d=_body(); s1=clean_seq(d.get('seq1','')); s2=clean_seq(d.get('seq2',''))
    if not s1 or not s2: return jsonify({'error':'Two sequences required'}),400
    if len(s1)>500 or len(s2)>500: return jsonify({'error':'Max 500 bp each'}),400
    st='PROTEIN' if detect_type(s1)=='PROTEIN' or detect_type(s2)=='PROTEIN' else 'DNA'
    fn=needleman_wunsch if d.get('method')=='global' else smith_waterman
    r=fn(s1,s2,st); r.update({'seq1_len':len(s1),'seq2_len':len(s2),'seq_type':st})
    return jsonify(r)

@app.route('/primers_advanced',methods=['POST'])
def primers_adv():
    d=_body(); raw=d.get('sequence','').strip()
    if not raw: return jsonify({'error':'No sequence'}),400
    seq=clean_seq(raw)
    if len(seq)<60: return jsonify({'error':'Min 60 bp required'}),400
    if len(seq)>50000: return jsonify({'error':'Max 50,000 bp for primer design'}),400
    return jsonify({'count':8,'pairs':advanced_primers(seq)})

_NT_DATABASES  = {'nt','refseq_rna','refseq_genomic'}
_PROT_DATABASES = {'nr','refseq_protein','swissprot','pdb'}

@app.route('/blast_submit',methods=['POST'])
def blast_sub():
    d=_body(); seq=clean_seq(d.get('sequence',''))
    if not seq: return jsonify({'error':'No sequence'}),400
    if len(seq)>10000: return jsonify({'error':'Max 10,000 bp for BLAST (NCBI limits apply)'}),400
    st=detect_type(seq); prog='blastp' if st=='PROTEIN' else 'blastn'
    # Default protein searches to refseq_protein (curated, ~20-60s) rather
    # than nr (exhaustive, 1-5+ min) when no database was specified at all —
    # nr is technically a *valid* protein database so the mismatch check
    # below wouldn't catch or correct this, and nr is explicitly the slow
    # option in our own UI, not what should be silently chosen by default.
    db=d.get('database','nt' if prog=='blastn' else 'refseq_protein')
    # The database dropdown holds both nucleotide and protein options
    # together, and nothing was resetting it when the pasted sequence
    # changed type — so a leftover protein selection (e.g. from a previous
    # "KRAS protein" sample) could get silently paired with a freshly
    # auto-detected DNA sequence. PROGRAM=blastn + DATABASE=refseq_protein
    # is not a valid combination for NCBI; it doesn't cleanly error, it
    # just behaves strangely (consistent with the wildly inflated RTOE seen
    # in practice). Enforce the pairing server-side regardless of what the
    # client sent, since this endpoint is the one place that actually knows
    # the true detected type.
    if prog=='blastn' and db not in _NT_DATABASES:
        db='nt'
    elif prog=='blastp' and db not in _PROT_DATABASES:
        db='refseq_protein'
    try: return jsonify(blast_submit(seq,program=prog,database=db))
    except Exception as e: return jsonify({'error':str(e)}),500

@app.route('/blast_poll',methods=['POST'])
def blast_p():
    rid=_body().get('rid','')
    if not rid: return jsonify({'error':'No RID'}),400
    try: return jsonify(blast_poll(rid))
    except Exception as e: return jsonify({'error':str(e)}),500

@app.route('/report',methods=['POST'])
def report():
    d=_body()
    if not d: return 'No data',400
    return Response(make_report(d),mimetype='text/html')

@app.route('/structure_search',methods=['POST'])
def structure_route():
    d=_body()
    gene=d.get('gene','').strip().upper()
    if not gene: return jsonify({'error':'No gene name'}),400
    if len(gene)>40 or not re.match(r'^[A-Z0-9\-\.]+$', gene):
        return jsonify({'error':'Invalid gene name format'}),400
    try: return jsonify(structure_search(gene))
    except Exception as e: return jsonify({'error':str(e)}),500

@app.route('/phylo_tree',methods=['POST'])
def phylo_tree_route():
    d=_body()
    raw=d.get('fasta','').strip()
    if not raw: return jsonify({'error':'No FASTA sequences provided'}),400
    seqs_list=parse_fasta(raw)
    if len(seqs_list)<2: return jsonify({'error':'Need at least 2 sequences in FASTA format'}),400
    # parse_fasta returns [{name, seq}] — build_phylo_tree needs {name: seq}
    seqs_dict={e['name']:e['seq'] for e in seqs_list}
    try: return jsonify(build_phylo_tree(seqs_dict))
    except Exception as e: return jsonify({'error':str(e)}),500

@app.route('/ramachandran',methods=['POST'])
def ramachandran_route():
    d=_body()
    pdb_url=d.get('pdb_url','').strip()
    if not pdb_url: return jsonify({'error':'No PDB URL provided'}),400
    try: return jsonify({'points': ramachandran_from_pdb(pdb_url)})
    except Exception as e: return jsonify({'error':str(e)}),500

@app.after_request
def add_headers(resp):
    resp.headers['X-Content-Type-Options'] = 'nosniff'
    resp.headers['X-Frame-Options'] = 'SAMEORIGIN'
    resp.headers['Referrer-Policy'] = 'strict-origin-when-cross-origin'
    resp.headers['Permissions-Policy'] = 'geolocation=(), microphone=(), camera=()'
    # Content-Security-Policy: restrict script/connection sources to exactly the
    # domains this app actually loads from. script-src blocks any injected
    # <script src="attacker.com/x.js"> from executing even if XSS gets a tag
    # onto the page. connect-src does the same for fetch()/XHR — matches the
    # server-side SSRF whitelist (_ALLOWED_FETCH_DOMAINS) as closely as possible,
    # plus the client-side PDB downloads the 3D viewer makes directly from the
    # browser (AlphaFold/RCSB), which never pass through our backend.
    resp.headers['Content-Security-Policy'] = (
        "default-src 'self'; "
        "script-src 'self' 'unsafe-inline' https://cdn.jsdelivr.net; "
        "style-src 'self' 'unsafe-inline' https://fonts.googleapis.com; "
        "img-src 'self' data:; "
        "font-src 'self' data: https://fonts.gstatic.com; "
        "connect-src 'self' https://*.ebi.ac.uk https://files.rcsb.org https://search.rcsb.org "
        "https://rest.uniprot.org https://blast.ncbi.nlm.nih.gov; "
        "frame-ancestors 'self'; "
        "object-src 'none'; base-uri 'self'"
    )
    # Gzip compress JSON/HTML/text responses over 1KB when the client supports it
    accept_enc = request.headers.get('Accept-Encoding', '')
    if ('gzip' in accept_enc and not resp.direct_passthrough
            and resp.content_type and
            any(t in resp.content_type for t in ('json', 'html', 'text'))
            and len(resp.get_data()) > 1024
            and 'Content-Encoding' not in resp.headers):
        buf = io.BytesIO()
        with gzip.GzipFile(fileobj=buf, mode='wb', compresslevel=6) as gz:
            gz.write(resp.get_data())
        resp.set_data(buf.getvalue())
        resp.headers['Content-Encoding'] = 'gzip'
        resp.headers['Content-Length'] = str(len(resp.get_data()))
        resp.headers['Vary'] = 'Accept-Encoding'
    return resp

if __name__=='__main__':
    print('\n'+'═'*52)
    print("  \U0001f9ec  Aaron's Archive - DNA Analysis Engine")
    print('═'*52)
    print('  Open: http://localhost:4040')
    print('  Stop: Ctrl+C')
    print('═'*52+'\n')
    app.run(debug=False,port=4040)
