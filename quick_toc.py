#!/usr/bin/python

import glob
import os

ignore_list = ("RAL_manual",)

def get_all_tex_files():
  res = []
  for root, dirnames, pathnames in os.walk('.'):
    to_add = glob.glob(os.path.join(root, "*.tex"))
    for file_name in to_add:
      if not any([f in file_name for f in ignore_list]):
        res.append(file_name)
  return res
  
def get_toc_formatted_line(line):
  def get_section_name(line):
    junk, name = line.split("{")
    return name.split("}")[0] # get rid of the rest of the line
  if '\\part' in line:
    return get_section_name(line)
  elif '\\chapter' in line:
    return '    '+get_section_name(line)
  elif '\\section' in line:
    return '        '+get_section_name(line)
  elif '\\subsection' in line:
    return '            '+get_section_name(line)
  elif '\\subsubsection' in line:
    return '                '+get_section_name(line)
  else:
    return None

def make_toc(file_names):
  res = []
  for path in file_names:
    with open(path, "r") as in_file:
      for line in in_file:
        formatted_line = get_toc_formatted_line(line)
        if formatted_line: res.append(formatted_line)
  return res

def main():
  file_names = get_all_tex_files()
  toc = make_toc(file_names)
  
  print "\n".join(toc)

if __name__=="__main__":
  main()


