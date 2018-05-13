import os, json, re, sqlite3, shutil, tempfile, math
# from .utils import createHTML, unique, open_file
from utils import createHTML, unique, open_file


def load_bed(filename):
  '''load a bed file into a list of lists'''
  con = open_file(filename)
  bed = list()
  for line in con:
    try: line = line.decode('utf-8')
    except: pass
    line = line.rstrip().split("\t")
    aux = []
    for i in range(len(line)):
      val = line[i]
      if i in [1,2,6,7,9]:
        val = int(val)
      aux.append(val)
    bed.append(aux)
  con.close()
  return bed


def get_assembly_from_fasta(fasta):
  '''return an assembly from a fasta file'''
  c = ''
  assembly = []
  for i in fasta:
    con = open_file(i)
    for line in con:
      try: line = line.decode('utf-8')
      except: pass
      line = line.strip()
      if line.startswith('>'):
        c = re.split(" |\|",re.sub("^>","",line))[0]
        assembly.append([c,0,0])
      else:
        assembly[-1][2] = assembly[-1][2] + len(line)
    con.close()
  return assembly


def get_assembly(assembly):
  '''return the specified assembly'''
  if assembly is None:
    return None
  if not isinstance(assembly,str):
    fasta = True
    valid = True
    for i in assembly:
      if not isinstance(i,str):
        fasta = False
        if len(i) != 3:
          valid = False
        elif not (isinstance(i[0],str) and isinstance(i[1],int) and isinstance(i[2],int)):
          valid = False 
      elif not os.path.isfile(i):
        fasta = False
    if fasta:
      return get_assembly_from_fasta(assembly)
    elif valid:
      return assembly
  if os.path.isfile(assembly):
    return load_bed(assembly)
  if assembly in ['NCBI36', 'GRCh37', 'GRCh38', 'GRCh37.bands', 'GRCh38.bands']:
    assembly = os.path.join(os.path.dirname(__file__), 'assemblies', assembly+'.json')
    return json.loads(open(assembly).read())
  print("Available assemblies are: 'NCBI36', 'GRCh37', 'GRCh38', 'GRCh37.bands', 'GRCh38.bands', but not '%s'." % assembly)
  return None


def chromosomesJSON(assembly):
  '''create json for chromosomes' cytobands'''
  d = dict()
  for i in range(len(assembly)):
    chromosome = assembly[i][0].lower()
    if chromosome not in d.keys():
      d[chromosome] = []
    d[chromosome].append(assembly[i][1:])

  return json.dumps(d)


def segmentation(track, cell):
  '''segmentation function'''
  seg = dict()
  for i in range(len(track)):
    try:
      value = float(track[i][4])
    except:
      continue
    j = (track[i][0],math.ceil(((track[i][1]+track[i][2])/2)/cell))
    try:
      seg[j][0] += 1
      seg[j][1] += value
    except:
      seg[j] = [1,value]

  track = []
  for i in seg:
    value = seg[i][1]/seg[i][0]
    track.append([i[0],(i[1]-1)*cell,i[1]*cell,None,value])

  return track


def genomemapJSON(data, assembly):
  '''create json for genome map'''
  d = dict()
  d['chromosomes'] = unique(map(lambda x: x[0], assembly))
  d['data'] = dict()
  d['max'] = max(map(lambda x: x[2], assembly))
  if not data is None:
    data = load_bed(data)
    if len(data)>100000:
      chr_max_len = dict()
      for i in d['chromosomes']:
        chr_max_len[i] = 0
      for i in range(len(assembly)):
        if assembly[i][2] > chr_max_len[assembly[i][0]]:
          chr_max_len[assembly[i][0]] = assembly[i][2]
      assembly_len = sum(chr_max_len.values())
      cell = math.ceil(assembly_len/100000)
      if cell>1:
        data = segmentation(data,cell)
    dataMin = float("inf")
    dataMax = float("-inf")
    for i in range(len(data)):
      if not data[i][0] in d['chromosomes']:
        continue
      value = None
      try:
        value = float(data[i][4])
      except:
        continue
      if data[i][0] not in d['data'].keys():
        d['data'][data[i][0]] = []
      d['data'][data[i][0]].append([data[i][1],data[i][2],value])
      if value < dataMin:
        dataMin = value
      if value > dataMax:
        dataMax = value
    d['dataDomain'] = [dataMin,dataMax]
  return json.dumps(d)


def genomemap(assembly, mapTrack = None, directory = "GenomeMap"):
  '''Create an interative genome map.
  
  Arguments:
      assembly -- a genome assembly, D3GB provides human assemblies ('NCBI36', 'GRCh37', 'GRCh38'), human assemblies with cytobands ('GRCh37.bands', 'GRCh38.bands'), or it can be an iterate of fasta file(s) paths or a single bed file path.
      mapTrack -- a bed file path with values to represent on the genome map. (default None)
      directory -- a string representing the directory where the graph will be saved. (default "GenomeBrowser")
  '''
  assembly = get_assembly(assembly)
  if assembly is not None:
    createHTML(directory, ["d3.min.js","jspdf.min.js","functions.js","genomemap.js"], genomemapJSON(mapTrack, assembly), chromosomesJSON(assembly))


def openDB(directory):
  '''open the sqlite db in the gb's directory'''
  return sqlite3.connect(os.path.join(directory, "Tracks.db"))


def createDB(directory):
  '''create sqlite to store tracks'''
  db = openDB(directory)
  c = db.cursor()
  c.execute("CREATE TABLE IF NOT EXISTS tbl_tracks (trackid INTEGER PRIMARY KEY, trackname TEXT, type TEXT, color TEXT, data TEXT)")
  c.execute("CREATE TABLE IF NOT EXISTS tbl_segments (trackid INTEGER, chr TEXT COLLATE NOCASE, start INTEGER, end INTEGER, name TEXT, score TEXT, strand TEXT, thickStart INTEGER, thickEnd INTEGER, itemRGB TEXT, blockCount INTEGER, blockSizes TEXT, blockStarts TEXT)")
  c.execute("CREATE INDEX location ON tbl_segments (chr,start,end)")
  db.commit()
  db.close()


def aux_genes(c):
  '''create the gene auxiliar database table'''
  c.execute("DROP TABLE IF EXISTS aux_genes")
  c.execute("CREATE TABLE aux_genes AS SELECT chr, start, end, name, score FROM tbl_segments NATURAL JOIN tbl_tracks WHERE (type='gene' OR type='exons')")


def add2DB(db, track, trackname, tracktype, color, scale):
  '''add tracks to a database connection'''
  if tracktype is None:
    tracktype = 'gene'
  if color is None:
    color = '#000'
  if tracktype not in ['value', 'score'] or scale is None or len(scale) != 2:
    scale = []
  try:
    scale = json.dumps(scale)
  except:
    scale = json.dumps([])

  c = db.cursor()
  c.execute('INSERT INTO tbl_tracks VALUES (NULL,?,?,?,?)',(trackname,tracktype,color,scale))
  trackid = c.lastrowid

  con = open_file(track)
  for line in con:
    try: line = line.decode('utf-8')
    except: pass
    line = line.rstrip().split("\t")
    row = [trackid]
    for i in range(12):
      val = None
      if i < len(line):
        if i in [1,2,6,7,9]:
          val = int(line[i])
        elif i == 4 and tracktype in ['value', 'score']:
          val = float(line[i])
        elif i == 5:
          if line[i] in ['+','-']:
            val = line[i]
        else:
          val = line[i]
      row.append(val)
    c.execute('INSERT INTO tbl_segments VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)',row)

  if tracktype in ['value', 'score'] and scale == '[]':
    c.execute("UPDATE tbl_tracks SET data=(SELECT '[' || min(cast(score as REAL)) || ',' || max(cast(score as REAL)) || ']' FROM tbl_segments WHERE trackid=?) WHERE trackid=?",(trackid,trackid))

  if tracktype in ['gene','exons']:
    aux_genes(c)

  db.commit()


def insert_track(c,uniq,track):
  '''check if track is in db and insert if not'''
  if not track in uniq.keys():
    tracktype = "gene" if track in ["CDS","gene"] else "domain"
    color = "cadetblue" if track in ["CDS","gene"] else "#000"
    c.execute("INSERT INTO tbl_tracks VALUES (NULL,?,?,?,?)",(track,tracktype,color,None))
    uniq[track] = c.lastrowid


class genomebrowser:

  def __init__(self, assembly, mapTrack = None, server = False, directory = 'GenomeBrowser'):
    '''Generates an interactive genome browser.
  
    Arguments:
      assembly -- a genome assembly, D3GB provides human assemblies ('NCBI36', 'GRCh37', 'GRCh38'), human assemblies with cytobands ('GRCh37.bands', 'GRCh38.bands'), or it can be an iterate of fasta file(s) paths or a single bed file path.
      mapTrack -- a bed file path with values to represent on the genome map. (default None)
      server -- a logical value to enable or disable the server mode. Server mode: designed to be shared as a website. Resulting folder should be added to the Apache applications directory and enable writting permissions. In this way the genome browser will be working as a web site. Local mode: the genome browser will be functional on your local machine. (default False)
      directory -- a string representing the directory where the graph will be saved. (default "GenomeBrowser")
    '''
    self.__directory__ = directory
    assembly = get_assembly(assembly)
    if assembly is not None:
      data = genomemapJSON(mapTrack, assembly)
      chromosomes = chromosomesJSON(assembly)
      if server:
        createHTML(directory, ["d3.min.js","jspdf.min.js","functions.js","images.js","genomebrowser.js","genomemap.js"], data, chromosomes)
        shutil.copy(os.path.join(os.path.dirname(__file__),'www','query.php'), directory)
      else:
        createHTML(directory, ["d3.min.js","jspdf.min.js","sql.js","functions.js","images.js","query.js","genomebrowser.js","genomemap.js"], data, chromosomes)
        print('Open the "index.html" file with Mozilla Firefox to see the graph. If you want to see this local file with other web browser, please visit the help section on the D3gb Web site http://d3gb.usal.es/')
      createDB(directory)
    else:
      self.__directory__ = None
      print('Invalid assembly.')


  def remove(self):
    '''Remove this genome browser'''
    shutil.rmtree(self.__directory__)


  def addTrack(self, track, trackname = None, tracktype = "gene", color = "#000", scale = None):
    '''Add tracks (bed files) to genome browser.
    
    Arguments:
      track -- a string representing the input bed file to be represented in the genome browser.
      trackname -- a string giving a name for the track.
      tracktype -- a string with the type of track should be drawn. Possible types are: "gene", "exons", "domain", "value" or "score". (default "gene")
      color -- a string giving the color of the track. (default "#000")
      scale -- a list with two values which specifies the minimun and maximun limits in the representation of the "score" or "value" tracks. By default, maximun and minimun scores are taken as the limits. (default None)
    '''
    if self.__directory__ is not None:
      if not tracktype in ['gene','domain','exons','value','score']:
        tracktype = 'gene'
      if not isinstance(trackname,str):
        trackname = os.path.split(track)[-1]
      db = openDB(self.__directory__)
      add2DB(db, track, trackname, tracktype, color, scale)
      db.close()


  def removeTrack(self,trackname):
    '''Removes a track from this genome browser by track name.
    
    Arguments:
      trackname -- a string giving the name of the track to remove.
    '''
    db = openDB(self.__directory__)
    c = db.cursor()
    ci = db.cursor()
    for row in ci.execute('SELECT trackid FROM tbl_tracks WHERE trackname=?',(trackname,)):      
      c.execute('DELETE FROM tbl_segments WHERE trackid=?',(row[0],))
      c.execute('DELETE FROM tbl_tracks WHERE trackid=?',(row[0],))
    aux_genes(c)
    db.commit()
    db.close()


  def addSequence(self, fastafile):
    '''Add sequences on Fasta format to genome browser.
    
    Arguments:
      fastafile -- a string or iterable representing the input Fasta file(s) to be added in the genome browser.
    '''
    directory = self.__directory__
    try: os.mkdir(os.path.join(directory,'sequences'))
    except: pass
    if isinstance(fastafile,str):
      fastafile = (fastafile,)
    for i in fastafile:
      con = open_file(i)
      seq = ''
      for line in con:
        try: line = line.decode('utf-8')
        except: pass
        line = line.rstrip()
        if line.startswith('>'):
          if not isinstance(seq,str):
            seq.close()
          chrom = re.split(' |\|',re.sub('^>','',line))[0]
          seq = open(os.path.join(directory,'sequences',chrom+'.fa'),'w')
          seq.write(">"+chrom+"\n")
        else:
          seq.write(line)
      seq.close()
      con.close()


  def addVCF(self, vcffile, trackname=None, show=None):
    '''Add vcf tracks to genome browser.
    
    Arguments:
      vcffile -- a string representing the input vcf file to be represented in the genome browser.
      trackname -- a string giving a name for the track.
      show -- an iterable giving the info features to display. (default None)
    '''

    def insert_track(c,name,master=False,json=None):
      tracktype = 'vcf' if master else 'vcfsample'
      c.execute('INSERT INTO tbl_tracks VALUES (NULL,?,?,?,?)',(name,tracktype,None,json))
      return c.lastrowid

    def insert_segment(c,trackid,chrom,pos,ID,info):
      c.execute('INSERT INTO tbl_segments (trackid,chr,start,end,name,score) VALUES (?,?,?,?,?,?)',(trackid,chrom,pos-1,pos,ID,info))

    def get_description(line):
      desc = []
      line = re.findall('<(.+)>',line)[0].split(',')
      for l in line:
        l = l.split('=')
        if l[0] == 'ID':
          desc.append(l[1])
        if l[0] == 'Description':
          desc.append(l[1].replace('"',''))
      return desc

    uniq_tracks = []
    trackid = 0
    chrom = ''
    pos = 0
    ID = ''
    info = ''
    data = {'info':{}}
    if show:
      data['show'] = show
    datacsq = []
    dataformat = {}
    db = openDB(self.__directory__)
    c = db.cursor()
    con = open_file(vcffile)
    if not isinstance(trackname,str):
      trackname = os.path.split(vcffile)[-1]

    for line in con:
      try: line = line.decode('utf-8')
      except: pass
      line  = line.rstrip()
      if line.startswith("##INFO=<ID="):
        desc = get_description(line)
        if desc[0] == 'CSQ':
          datacsq = re.findall('Format: (.+)$',desc[1])[0].split('|')
        else:
          data['info'][desc[0]] = desc[1]
      elif line.startswith("##FORMAT=<ID="):
        desc = get_description(line)
        dataformat[desc[0]] = desc[1]
      elif line.startswith("#CHROM"):
        data = {} if len(data['info']) is 0 else data
        uniq_tracks.append(insert_track(c,trackname,True,json.dumps(data)))
        aux = line.split("\t")
        for i in range(9,len(aux)):
          uniq_tracks.append(insert_track(c,aux[i],False,json.dumps(dataformat)))
      elif not line.startswith("#"):
        aux = line.split("\t")
        trackid = uniq_tracks[0]
        chrom = aux[0]
        pos = int(aux[1])
        ID = aux[2] if aux[2]!='.' else None
        score = ["REF="+aux[3]]
        if aux[4]!='.':
          score.append("ALT="+aux[4]);
        if aux[5]!='.':
          score.append("QUAL="+aux[5]);
        try:
          info = aux[7].split(';')
          for infi in info:
            if infi.startswith("CSQ="):
              csq = infi.replace('CSQ=','').split(',')[0].split('|')
              for i in range(len(datacsq)):
                if csq[i]!="":
                  score.append(datacsq[i]+"="+csq[i])
            else:
              score.append(infi)
        except:
          pass
        info = '|'.join(score)

        insert_segment(c,trackid,chrom,pos,ID,info)
        try: fmt = aux[8].split(":")
        except: continue
        for i in range(1,len(uniq_tracks)):
          attr = aux[i+8].split(":")
          if len(fmt) == len(attr):
            trackid = uniq_tracks[i]
            for j in range(len(fmt)):
              attr[j] = fmt[j]+'='+attr[j]
            info = '|'.join(attr)
            insert_segment(c,trackid,chrom,pos,ID,info)

    con.close()
    db.commit()
    db.close()


  def addGFF(self,gfffile):
    '''Add tracks in a gff file to  genome browser.
    
    Arguments:
      gfffile -- a string representing the input gff file to be represented in the genome browser.
    '''

    uniq_tracks = dict()
    trackName = ''
    IDs = dict()
    segments = []

    scaffold = ''
    start = 0
    end = 0
    name = ''
    score = ''
    strand = None
    thickStart = None
    thickEnd = None

    db = openDB(self.__directory__)
    c = db.cursor()
    con = open_file(gfffile)
    for line in con:
      try: line = line.decode('utf-8')
      except: pass
      line  = line.rstrip()
      if line.startswith('#') or line is '':
        continue
      aux = line.split("\t")
      trackName = aux[2]
      try:
        regex = re.compile(r"\b(\w+)\=([^\=]*)(?=;\w+\=|$)")
        attrs =  dict(regex.findall(aux[8]))
      except:
        attrs = dict()
      start = int(aux[3])-1
      end = int(aux[4])
      try: name = attrs.pop('Name')
      except: name = None
      try: ID = attrs.pop('ID')
      except: ID = None
      if (name is None) and not (ID is None):
        name = ID
      if trackName == 'exon':
        try: name = attrs.pop('Parent')
        except: continue
        size = str(end - start)
        # get start, blockCount, blockSizes, blockStarts from parent:
        f = [segments[IDs[name]][x] for x in [0,2,9,10,11]] 
        st = str(start - f[1])
        segments[IDs[name]][9] = f[2] + 1 # blockCount
        segments[IDs[name]][10] = f[3] + size + ',' # blockSizes
        segments[IDs[name]][11] = f[4] + st + ',' # blockStarts
        # change track type to exons:
        if segments[IDs[name]][9] == 1: 
          c.execute('UPDATE tbl_tracks SET type=?,color=? WHERE trackid=?',('exons','goldenrod',f[0]))
        continue
      insert_track(c,uniq_tracks,trackName)
      score = json.dumps(attrs,separators=('|','=')) if len(attrs)!=0 else None
      score = re.sub('["{}]','',score)
      scaffold = aux[0]
      strand = aux[6] if aux[6] in ['+','-'] else None
      thickStart = None
      thickEnd = None
      if aux[7] in ['0','1','2']:
        if strand!='-':
          thickStart = start+int(aux[7])
          thickEnd = end
        else:
          thickStart = start
          thickEnd = end-int(aux[7])
      segments.append([uniq_tracks[trackName],scaffold,start,end,name,score,strand,thickStart,thickEnd,0,"",""])
      if not ID is None and not ID in IDs.keys():
        IDs[ID] = len(segments)-1

    con.close()
    cur = db.cursor()
    for row in segments:
      c.execute('INSERT INTO tbl_segments (trackid,chr,start,end,name,score,strand,thickStart,thickEnd,blockCount,blockSizes,blockStarts) VALUES (?,?,?,?,?,?,?,?,?,?,?,?)',row)
      
    aux_genes(c)
    db.commit()
    db.close()


def gbk2genomebrowser(gbkfile, server = False, directory = "GenomeBrowser"):
  '''Creates an interactive genome browser from a GenBank file.
  
  Arguments:
    gbkfile -- a string representing the input GenBank file to be represented in the genome browser.
    server -- a logical value to enable or disable the server mode. Server mode: designed to be shared as a website. Resulting folder should be added to the Apache applications directory and enable writting permissions. In this way the genome browser will be working as a web site. Local mode: the genome browser will be functional on your local machine. (default False)
    directory -- a string representing the directory where the graph will be saved. (default "GenomeBrowser")
  '''
  current = 1
  track = ''
  string = ''
  sequence = ''

  uniq_tracks = dict()
  assembly = []

  tmp = tempfile.mkdtemp()
  createDB(tmp)
  db = openDB(tmp)
  c = db.cursor()

  def trackData(c,string):
    if string != '':
      string = string.split('\n/')
      scaffold = assembly[-1][0]
      start = 0
      end = 0
      name = None
      score = None
      strand = '+'
      blockCount = None
      blockSizes = None
      blockStarts = None
      pos = string.pop(0)
      pos = re.sub('<|>|\\n','',pos)
      if 'complement' in pos:
        strand = '-'
        pos = re.sub('complement\(|\)','',pos)
      if 'join' in pos:
        c.execute('UPDATE tbl_tracks SET type=?,color=? WHERE trackid=?',('exons','goldenrod',uniq_tracks[track]))
        pos = re.sub('join\(|\)','',pos)
        pos = pos.split(",")
        blockCount = len(pos)
        pos = list(map(lambda x: list(map(int,x.split('..'))),pos))
        pos.sort()
        start = pos[0][0]-1
        end = pos[-1][1]
        blockSizes = ','.join(list(map(lambda x: str(x[1]-(x[0]-1)), pos)))
        blockStarts = ','.join(list(map(lambda x: str((x[0]-1)-start), pos)))
      else:
        pos = list(map(int,pos.split("..")))
        start = pos[0]-1
        end = pos[1]
      for s in string:
        if s.startswith("gene="):
          name = re.sub('gene=|"','',s)
        elif s.startswith("locus_tag=") and name is None:
          name = re.sub('locus_tag=|"','',s)
        elif not s.startswith("translation="):
          if score is None:
            score = s.replace('"','')
          else:
            score = score + '|' + s.replace('"','')
      c.execute('INSERT INTO tbl_segments (trackid,chr,start,end,name,score,strand,blockCount,blockSizes,blockStarts) VALUES (?,?,?,?,?,?,?,?,?,?)',(uniq_tracks[track],scaffold,start,end,name,score,strand,blockCount,blockSizes,blockStarts))

  con = open_file(gbkfile)
  for line in con:
    try: line = line.decode('utf-8')
    except: pass
    aux = re.split(' +',line.strip())
    if aux[0]=="LOCUS":
      assembly.append([aux[1],0,int(aux[2])])
      current = 1
    elif aux[0]=="VERSION":
      if aux[1].find(assembly[-1][0]) == 0:
        assembly[-1][0] = aux[1]
      current = 1
    elif aux[0]=="FEATURES":
      current = 2
    elif aux[0]=="ORIGIN":
      trackData(c,string)
      string = ''
      current = 3
      sequence = open(os.path.join(tmp,assembly[-1][0]+".fa"),'w')
      sequence.write(">"+assembly[-1][0]+"\n")
    elif re.search("^[A-Z]",line[0]):
      current = 1
    else:
      if current==2:
        if line[5]!=" ":
          trackData(c,string)
          string = ''
          track = aux[0]
          if(track!="source"):
            string = aux[1]
            insert_track(c,uniq_tracks,track)
        else:
          if(track!="source"):
            string = string+'\n'+line.strip()
      if current==3:
        if aux[0]!="//":
          sequence.write(''.join(aux[1:]))
        else:
          sequence.close()
  con.close()
  aux_genes(c)
  db.commit()
  db.close()

  gb = genomebrowser(assembly,None,server,directory)
  os.mkdir(os.path.join(directory,"sequences"))
  for chrom in assembly:
    shutil.move(os.path.join(tmp,chrom[0]+".fa"),os.path.join(directory,"sequences",chrom[0]+".fa"))
  os.remove(os.path.join(directory,'Tracks.db'))
  shutil.move(os.path.join(tmp,'Tracks.db'),os.path.join(directory,'Tracks.db'))
  shutil.rmtree(tmp)
  return gb
