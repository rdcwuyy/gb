import os, shutil, re, webbrowser, gzip

def open_file(filename):
  con = open(filename,'rb')
  if con.read(1) == b'\x1f' and con.read(1) == b'\x8b':
    return gzip.open(filename,'rb')
  else:
    return open(filename,'r')

def unique(items):
  uniq = list()
  for i in items:
    if not i in uniq:
      uniq.append(i)
  return uniq

www = os.path.join(os.path.dirname(__file__), 'www')

def createHTML(directory, dependencies, data, chromosomes):
  if os.path.exists(directory):
    shutil.rmtree(directory)
  os.makedirs(directory)
  html = open(os.path.join(www, 'template.html')).read()
  name = directory.split(os.sep)[-1]
  html = html.replace('<!--title-->', name)
  dep = '<!--head-->'
  for depend in dependencies:
    dirName = ''
    if re.search(".css$",depend):
      dep = dep + "\n" + '<link rel="stylesheet" type="text/css" href="styles/' + depend + '"></link>'
      dirName = os.path.join(directory,'styles')
    else:
      dep = dep + "\n" + '<script type="text/javascript" src="scripts/' + depend + '"></script>'
      dirName = os.path.join(directory,'scripts')
    if not os.path.exists(dirName):
      os.mkdir(dirName)
    shutil.copy(os.path.join(www,depend), dirName)

  html = html.replace('<!--head-->', dep)
  chromosomes = '<script type="application/json" id="chromosomes">' + chromosomes + '</script>'
  html = html.replace('<!--body-->','<!--body-->\n' + chromosomes)
  data = '<script type="application/json" id="data">' + data + '</script>'
  html = html.replace('<!--body-->','<!--body-->\n' + data)
  index_path = os.path.join(directory, "index.html")
  con = open(index_path,"w")
  con.write(html)
  con.close()
  print("The graph has been generated in the '%s' folder." % directory)
