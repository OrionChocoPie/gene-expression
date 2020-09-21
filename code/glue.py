from __future__ import print_function
from subprocess import call
from math import ceil
import os
import shutil

def label():
  global filenames;
  newnames = [];
  i = 0;
  for filename in filenames:
    newFilename = "{}/{}.png".format(temppath,alphabet[i]);
    newnames.append(newFilename);
    argv = ["convert",filename,"-fill","black","-pointsize","200","-gravity","northwest","-annotate","+100+100",alphabet[i],newFilename];
    print("labeling:",filename);
    print (argv);
    call(argv);
    i += 1;
  filenames = newnames;

filenames = [
  "./QN_01_rect.png",
  "./DESeq2_01_rect.png",
  "./Sh00_01_rect.png",
  "./Sh01_01_rect.png",
  "./Sh10_01_rect.png",
  "./Sh11_01_rect.png",
  "./Sh11_01_rect.png",
  "./Sh_571_0_01_rect.png",
  "./Sh_571_1_01_rect.png",
  "./Sh_10558_0_01_rect.png",
  "./Sh_10558_1_01_rect.png",
  "./Sh_11154_0_01_rect.png",
  "./Sh_11154_1_01_rect.png",
  "./Sh_AGL_0_01_rect.png",
  "./Sh_AGL_1_01_rect.png",
  "./Sh_75_0_01_rect.png",
  "./Sh_75_1_01_rect.png",
  "./legendMM.png",
  "./vio16SL.png",
  ];
targetFilename = "./all16.png";
cols = 4;

rows = int(ceil(len(filenames)/cols));
temppath = "./temp";
alphabet = "abcdefghijklmnopqrstuvwxy0123456789";
tempfiles = [];

if not os.path.exists(temppath):
  os.makedirs(temppath);
else:
  print("temppath exists");

label();#if needed

for k in range(0,rows):
  offset = k*cols;
  tempfile = "{}/_{}.png".format(temppath,str(k));  
  tempfiles.append(tempfile);
  names = filenames[offset:offset+cols]
  argv = ["convert"]+names+["+append"]+[tempfile];
  print (argv);
  print("gluing:"," ".join(names));
  call(argv);

argv = ["convert"]+tempfiles+["-append"]+[targetFilename];

print("gluing:"," ".join(tempfiles));
print (argv);
call(argv);

shutil.rmtree(temppath)
