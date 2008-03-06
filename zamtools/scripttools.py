#!/usr/bin/env python

#LAST MODIFIED: 07-25-06

import sys, os, glob


def ExpandArgs(Args):
  """Expands arguments that contains wildcards.
On Linux systems, does nothing since command interpreter
already does this."""
  if os.name == "nt":
    NewArgs = []
    for a in Args:
      if ("*" in a or "?" in a) and not a.startswith("--"):
        NewArgs.extend(glob.glob(a))
      else:
        NewArgs.append(a)
    return NewArgs
  else:
    return Args

def ArgAsType(Val, Template):
  """Returns Val as the type of Template"""
  #for none we will return True
  if Template is None:
    return True
  else:
    t = type(Template)
    return t(Val)
      

def ParseArgs(Argv, Defaults = {}):
  """Parses a list of arguments into a dictionary.
The keys 0,1,..n give the consecutive non-option arguments.
The key "FLAGS" gives a list of flags used.
The key "ARGS" gives a list of non-flag arguments.
The key "NARG" gives the number of arguments.
Other keys give the values of flag arguments.
Flags can be '-x', '--word', or '--word=value'."""
  ind = -1
  Ret = Defaults.copy()
  Ret["FLAGS"] = []
  Argv = ExpandArgs(Argv)
  for i in range(0, len(Argv)):
    a = Argv[i]
    if a.startswith("--"):
      #word flags
      if "=" in a:
        Flag = a[2:].split("=")[0]
        Val = a[2:].split("=")[1]
      else:
        #get the following arg
        Flag = a[2:]
        if i + 1 < len(Argv):
          Val = Argv[i+1]
        else:
          Val = ""
      #try to get the same type as defaults
      if Flag in Defaults:
        Ret[Flag] = ArgAsType(Val, Defaults[Flag])
      else:
        Ret[Flag] = Val
      Ret["FLAGS"].append(Flag)
    elif a.startswith("-") and len(a) > 1 and not a[1] in "0123456789":
      #one letter flags, use following arg as value
      if i + 1 < len(Argv):
        Val = Argv[i+1]
      else:
        Val = ""
      for Flag in a[1:]:
        #try to get the same type as defaults
        if Flag in Defaults:
          Ret[Flag] = ArgAsType(Val, Defaults[Flag])
        else:
          Ret[Flag] = Val
        Ret["FLAGS"].append(Flag)
    else:
      #normal argument
      ind += 1
      Ret[ind] = a
  Ret["NARG"] = ind + 1
  Ret["ARGS"] = [Ret[i] for i in range(0, ind+1)]
  return Ret
