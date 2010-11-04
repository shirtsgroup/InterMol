import os, string, pdb, copy, re
from collections import deque

class GromacsTopologyParser:
	def __init__(self, topfile=None, defines=None):
		"""A simple class to handle the parsing of Gromacs topology files
		
		Generates an instance of this object preprocessing topfile and adding
		the specified defines
		
		Args:
			topfile: The name of a topology file to parse
			defines: a set of default defines to add
		"""	
		self.includes = set()		# set storing includes
		self.expanded = list()		# list of expanded lines
		self.defines = dict()		# list of defines
		if defines:
			self.defines.union(defines)
		self.defines["FLEX_SPC"] = None
		self.defines["POSRE"] = None
		if topfile:
			self.preprocess(topfile)		
		self.comments = list()
		self.directives = list()	# list of directives				
		self.systems = list() 		# list of system directives
		self.molecules = list()		# list of molecule directives
		self.parameters = list()	# list of parameter directives
	
	def read(self, verbose = False):
		"""reads information from a previously preprocessed topology file into
		a gromacstopologyfile object
		
		Args:
			verbose: verbose output
		"""
		sysDirective = re.compile(r"""
		  \[[ ]{1}
		  ((?P<system>molecules|system)				# matches molecules or system 
		  |
		  (?P<parameter>defaults|.*types|nonbond_params)	# matches .*atomtypes
		  |
		  (?P<molecule>[\w]+))		 			# everything else alphanumeric
		  [ ]{1} \]
		""", re.VERBOSE)
		for line in self.expanded:
			match = sysDirective.match(line)
			if match:
				if verbose:
					print match.groups()
				directive = self.Directive(match.group(0), header='')
				if match.group('system'):
					self.systems.append(directive)
				elif match.group('parameter'):
					self.parameters.append(directive)
				elif match.group('molecule'):
					self.molecules.append(directive)	
			else:
				directive.lines.append(line)

	def preprocess(self, filename, verbose = False):
		"""preprocess a topology file
		
		Preprocesses for #include, #define, #ifdef directives and handles them
		
		Args:
			filename: the filename to preprocess
			verbose: verbose output
		"""
		lineReg = re.compile(r"""
		 [ ]*					# omit spaces
		 (?P<directive>[^\\;\n]*) 	# search for directive section
		 [ ]*					# omit spaces
		 (?P<ignore_newline>\\)?		# look for the \ for newline
		 (?P<comment>;.*)?			# look for a comment
		""", re.VERBOSE)
		preReg = re.compile(r"""
		 (?:\#include[ ]+)
		  (?:"|<)
		  (?P<include>.+\.itp|.+\.top)
		  (?:"|>)
	 	 |	
		 (?:\#define[ ]+)
		  (?P<defineLiteral>[\S]+)
		  (?P<defineData>[ ]+.+)?
		 |
		 (?:\#undef[ ]+)
		  (?P<undef>[\S]+)
		 |
		 (?:\#ifdef[ ]+)
		  (?P<ifdef>[\S]+)
		 |
		 (?:\#ifndef[ ]+)
		  (?P<ifndef>[\S]+)
		 |
		 (?P<else>\#else)
		 |
		 (?P<endif>\#endif)
		""", re.VERBOSE)
		condDepth = 0
		write = [True]	
		fd = self.open(filename)
		if fd != None:
			lines = deque(fd)	# initial read of root file
			while(lines):		# while the queue isn't empty continue preprocessing
				line = lines.popleft()
				if verbose:
					print "=========================" 
					print "Original line:", line
				match = lineReg.match(line)	# ensure that the line matches what we expect!
				if match:
					line = match.group('directive')
					comment = match.group('comment')
					if verbose:
						print "Directive:", line
						print "Comment:", comment
					while True:		# expande the lines
						if match.group('ignore_newline'):
							if verbose: print "Ignore newline found"
							if lines:
								temp = lines.popleft()
								match = lineReg.match(temp)
								line += match.group('directive')
							else:
								print "WARNING: Previous line continues yet EOF encountered"
								break	# EOF so we can't loop again		
						else:
							break	# line terminates	
					line = line.strip()	# just in case theres whitespace we missed
						
					match = preReg.match(line)
					if match != None:
						if match.group('ifdef') and write[condDepth]:
							if verbose:
								print 'Found an ifdef:', match.group('ifdef')
							condDepth += 1
							ifdef = match.group('ifdef')
							if ifdef in self.defines:
								write.append(True)
							else:
								write.append(False)	
						elif match.group('ifndef') and write[condDepth]:
							if verbose:
								print 'Found an ifndef:', match.group('ifndef')
						 	condDepth += 1
							ifndef = match.group('ifndef')
							if ifndef not in self.defines:
								write.append(True)
							else:
								write.append(False)	
						elif match.group('else'):
							if verbose:
								print 'Found an else'
							if condDepth > 0:
								write[condDepth] = not write[condDepth]
							else:
								print "WARNING: Found an else not associated with a conditional"
						elif match.group('endif'):
							if verbose:
								print 'Found an endif'
							if condDepth > 0:
								condDepth -= 1
								write.pop()
							else:
								print "ERROR: Found an endif not associated with a conditional"
						elif write[condDepth]:
							if match.group('include'):
								if verbose:
									print "Found a include:", line
								include = match.group('include')
								fd = self.open(include)
								if fd != None:
									tempList = list(fd)
									tempList.reverse()
									lines.extendleft(tempList)
							elif match.group('defineLiteral'):
								if verbose:
									print "Found a define:", line
								define = match.group('defineLiteral')
								if define not in self.defines:
									self.defines[define] = match.group('defineData')
								else:
									print "WARNING: Overriding define:", define
									self.define[define] = match.group('defineData')
							elif match.group('undef'):
								if verbose:
									print "Found a undefine:", line
								undef = match.group('undef')
								if undef in self.defines:
									self.defines.pop(undef)		
					elif write[condDepth]:
						if line != '':
							for define in self.defines:
								if define in line:
									line = line.replace(define, self.defines[define])
							if verbose:
								print "Writing:", line
							self.expanded.append(line)
						if comment:
							self.comments.append(comment)
				else:
					print "ERROR: Unreadable line!"	
	
	def open(self, filename):
		"""Open a file and add to includes list
	
		Checks if a file exists in the current directory or in the GMXLIB environment then either opens or fails
		
		Args:
			filename: the name of the file to open
	
		Returns:
			file descriptor of the file to be openend
		
		Raises:
			IOError when the file cannot be opened for any reason 
		"""
		if filename in self.includes:
			print "WARNING: Omitting file ", filename, ". It has already been included!"
			return None
		temp = filename
		if not os.path.exists(filename):
			os.path.join(os.environ['GMXLIB'], filename)
		try:
			fd = open(filename)
			self.includes.add(temp)
			return fd
		except IOError:
			print "WARNING: File ", filename, " could not be found!"
			return None
	
	def write(self, filename):
		fout = open(filename, 'w')
		for line in self.expanded:
			fout.write(line)
		fout.close()
		return

	def getComments(self):
		return self.comments
	
	def getExpanded(self):
		return self.expanded

	def getIncludes(self):
		return self.includes
	
	def getDefines(self):
		return self.defines
	
	#=============================================================================================
	# Containers for *.top file DIRECTIVES, e.g.: [ defaults ], [ atomtypes ]
	#=============================================================================================
	class Directive(object):
		"""A base class for for the Directive container"""

		def __init__(self, name, header=None):
       			self.classType = 'Directive'
       		 	self.name = name        # Can be 'name', or be '[ name ]\n' (It will be fixed to be the latter).
			self.header = header    # any "; <comment>\n" text between the name and data (can have multiple newlines).
        		self.lines = []
	
                def fixName(self):
                        """if name is 'name', convert it to '[ name ]' """
                        self.name = '[ ' + self.name.replace('[','').replace(']','').strip() + ' ]\n'
                        return

                def getName(self):
                        return self.name

                def setName(self, name):
                        self.name = name
                        return

                def getHeader(self):
                        return self.header

                def setHeader(self, header):
                        self.header = header
                        return

                def __repr__(self):
                        """Returns a string representation that can be printed or written to file."""

                        outtxt = self.name + '\n' + self.header
                        for line in self.lines:
                                outtxt += line
                        outtxt += '\n'  # each directive needs a trailing newline for coding
                        return outtxt
