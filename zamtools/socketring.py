#!/usr/bin/env python

#LAST MODIFIED: 09-10-06


import sys, os, time, copy, random, cPickle, shutil
import socket, select 


#LOCAL VARIABLES WHICH HAVE SPECIAL MEANING IN SOCK EXECS:
#SockData: data that was sent over with the exec command
#SockResult: data that is to be sent back after an exec
#SockCloseExec: string of commands to be executed when sock closes
#SockTempFiles: list of temporary files or paths to be removed upon close


#GLOBALS

#socket message constants
MsgLen = 4096
MsgEnd = "[[MSGEND]]"
MsgDone = "[[DONE]]"
MsgExec = "[[EXEC]]"
MsgEval = "[[EVAL]]"
MsgRelease = "[[RELEASE]]"
MsgFinish = "[[FINISH]]"

#random connect wait
RandTime = 5.
#maximum size for the file sockets.txt, in bytes
MaxStatsSize = 1024*1024
#maximum number of nodes that can connect
MaxConnectDflt = 256
#maximum time a client node can spend idle, in sec
MaxIdleDflt = 120.
#amount of time to spend between checking sockets
IntervalDflt = 5.

#port defaults
PortMinDflt = 9915
PortMaxDflt = 10015

#indicate when a task is run on a client in verbose mode
TaskSignal = False

#blocking event
BlockingEvent = 10035


#======== ACCESSORY FUNCTIONS ========

def GetPort(PortMin = PortMinDflt, PortMax = PortMaxDflt):
  """Returns the first available socket port in a range."""
  for PortNum in range(PortMin, PortMax):
    #check to see if port is bound
    try:
      sc = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
      sc.bind((socket.gethostname(), PortNum))
    except socket.error:
      Bound = True
    else:
      Bound = False
    sc.close()
    if not Bound:
      return PortNum
  return 0



#======== SOCKET WRAPPER ========

class SSockClass:

  def __init__(self, Sock = None):
    """A sockets wrapper class that manages simple sends and
    receives, with the ability to send exec or eval commands."""
    self.Buf = ""
    self.Busy = False
    self.Addr = ""
    self.Port = 0
    self.Name = ""
    self.Connected = False
    self.Closed = False
    self.Task = None
    self.TimeUp = None
    self.TimeStart = None
    self.TimeStop = None
    self.SecWork = 0.
    self.SecIdle = 0.
    #blank variable with commands to be executed upon finishing
    self.SockCloseExec = ""
    self.SockTempFiles = []
    if Sock is None:
      self.GetNewSock()
    else:
      self.Sock = Sock
      self.Connected = True
      self.__GetNames()  

  def __WorkTimeStart(self):
    "Indicates start of working on exec/eval."
    self.TimeStart = time.time()
    if not self.TimeStop is None:
      self.SecIdle += self.TimeStart - self.TimeStop

  def __WorkTimeStop(self):
    "Indicates stop of working on exec/eval."
    self.TimeStop = time.time()
    if not self.TimeStart is None:
      self.SecWork += self.TimeStop - self.TimeStart

  def ElapsedSec(self):
    "Returns the number of seconds this process has been running."
    if not self.TimeStart is None:
      return max(time.time() - self.TimeStart, 0)
    else:
      return 0

  def GetNewSock(self):
    """Creates a new socket."""
    self.Sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    self.Closed = False  

  def __GetNames(self):
    """Gets the socket properties."""
    self.Addr = ""
    self.Port = 0
    try:
      self.Addr, self.Port = self.Sock.getpeername()
    except:
      try:
        self.Addr, self.Port = self.Sock.getsockname()
      except:
        self.Addr = "[addr error]"  
    self.Name = ""
    try:
      self.Name = socket.gethostbyaddr(self.Addr)[0]
    except:
      self.Name = "[name error]"    

  def Listen(self, Port, MaxConnect = MaxConnectDflt):
    """Starts a listening socket."""
    try:
      self.Sock.setblocking(0)
      #bind the socket to us
      self.Sock.bind((socket.gethostname(), Port))
      self.__GetNames()
      #start listening on that port
      self.Sock.listen(MaxConnect)
    except socket.error, errmsg:
      self.Connected = False
      return False
    else:
      self.Connected = True
      return True
    
  def Connect(self, Host, Port):
    """Connects the socket to the Host on Port."""
    try:
      self.Sock.connect((Host, Port))
    except socket.error, errmsg:
      self.Connected = False
      return False
    else:
      self.TimeUp = time.time()
      self.Connected = True
      self.__GetNames()
      return True    

  def __DelTempFiles(self):
    "Removes any temporary files."
    for f in self.SockTempFiles:
      if os.path.isfile(f):
        os.remove(f)
      elif os.path.isdir(f):
        shutil.rmtree(f)
    self.SockTempFiles = []  

  def Close(self):
    "Closes the socket and cleans up."
    self.Sock.close()
    self.Connected = False
    self.Closed = True
    self.Busy = False
    self.TimeStart, self.TimeStop = None, None
    self.Buf = ""
    #run any optional tidying up tasks
    if len(self.SockCloseExec) > 0:
      exec(self.SockCloseExec)
      self.SockCloseExec = ""
    #check for temp files to delete
    self.__DelTempFiles()

  def Send(self, Msg):
    "Sends a string message through the socket."
    s = Msg + MsgEnd
    try:
      while len(s) > 0:
        n = self.Sock.send(s)
        s = s[n:]
    except socket.error, errmsg:
      self.Close()
      return False
    return True

  def Recv(self):
    "Receives a string message."
    #keep receiving until the full msg arrives
    LastInd = 0
    while not MsgEnd in self.Buf[LastInd:]:
      try:
        t = self.Sock.recv(MsgLen)
      #if there is an error, return None
      except socket.error, errmsg:
        #skip over blocking warning
        if not errmsg[0] == BlockingEvent:
          return None
      else:
        self.Buf += t
      if len(t) == 0:
        return None
      LastInd = max(len(self.Buf) - len(t) - len(MsgEnd), 0)
    l = LastInd + self.Buf[LastInd:].find(MsgEnd)
    s = self.Buf[0:l]
    self.Buf = self.Buf[l+len(MsgEnd):]
    #we need to check if this is a special message:
    if s.startswith(MsgExec):
      #get the lengths
      l1 = len(MsgExec)
      l2 = l1 + 10
      l3 = int(s[l1:l2]) + l2
      #this is a command to execute
      Cmd = s[l2:l3]
      #get the data
      SockData = cPickle.loads(s[l3:])
      #make the result variable
      SockResult = None
      #make a local socket closing variable so it can be modified
      SockCloseExec = self.SockCloseExec
      SockTempFiles = self.SockTempFiles
      #execute it
      exec(Cmd)
      #check to see if there was a result
      Ret = MsgDone
      if not SockResult is None:
        Ret += cPickle.dumps(SockResult)
      #check if the socket closing variable has been modified
      self.SockCloseExec = SockCloseExec
      self.SockTempFiles = SockTempFiles
      #send back a done string
      if self.Send(Ret):
        return MsgExec
      else:
        #error so close the socket
        self.Close()
        return None
    elif s.startswith(MsgEval):
      #get the lengths
      l1 = len(MsgEval)
      l2 = l1 + 10
      l3 = int(s[l1:l2]) + l2
      #this is a command to evaluate
      Cmd = s[l2:l3]
      #get the data
      SockData = cPickle.loads(s[l3:])
      #execute it
      Ret = MsgDone + cPickle.dumps(eval(Cmd))
      #send back a done string
      if self.Send(Ret):
        return MsgEval
      else:
        #error so close the socket
        self.Close()
        return None
    elif s.startswith(MsgDone):
      #this is a done message
      self.Busy = False
      #see if there is a return value
      Ret = s[len(MsgDone):]
      self.__WorkTimeStop()
      if len(Ret) > 0:
        return cPickle.loads(Ret)
      else:
        return MsgDone
    else:
      #this is a plain old message
      return s

  def SendExec(self, ExecStr, Data = None):
    "Sends a string command to be executed by the processor."
    self.Busy = True
    self.__WorkTimeStart()
    l = len(ExecStr)
    DataStr = cPickle.dumps(Data)
    return self.Send(MsgExec  + ("%010d" % l) + ExecStr + DataStr)

  def SendEval(self, EvalStr, Data = None):
    "Sends a string command to be evaluated by the processor."
    self.Busy = True
    self.__WorkTimeStart()
    l = len(EvalStr)
    DataStr = cPickle.dumps(Data)
    return self.Send(MsgEval + ("%010d" % l) + EvalStr + DataStr)

  def FracWork(self):
    "Returns the fraction of the time this processor is working."
    Tot = self.SecIdle + self.SecWork
    if Tot > 0.:
      return self.SecWork / Tot
    else:
      return 1.



#======== TASK CLASS ========

class TaskClass:
  def __init__(self, TaskStr, ID, Data = None, Timeout = None):
    self.TaskStr = TaskStr
    self.ID = ID
    self.Done = False
    self.Client = None
    self.Data = Data
    self.Result = None
    self.TimeStart = time.time()
    self.TimeStop = None
    self.Timeout = Timeout
  def ElapsedTime(self):
    if self.TimeStop is None or self.TimeStart is None:
      return None
    else:
      return self.TimeStop - self.TimeStart
  


#======== SOCKET RING CLASS ========

class RingClass:

  def __init__(self, Path = ".", Verbose = True, MaxConnect = MaxConnectDflt,
    PortMin = PortMinDflt, PortMax = PortMaxDflt, Interval = IntervalDflt,
    MaxIdle = MaxIdleDflt, MaxNeeded = 0, ClientTimeout = None):
    #list out the variables
    self.Path = Path
    self.MaxConnect = MaxConnect
    self.IsServer = False
    self.Server = ""
    self.Port = 0
    self.Verbose = Verbose
    self.Interval = Interval
    self.MaxIdle = MaxIdle
    self.LastID = 0
    self.Tasks = []
    self.MaxNeeded = MaxNeeded
    self.ClientTimeout = ClientTimeout
    
    #check to see if we need to make the path
    if not os.path.isdir(self.Path): os.mkdir(self.Path)
    
    #DETERMINE SERVER/CLIENT ROLES
    self.__DecideServer()
    if self.Server == "":
      print "Problem determining main server."
    elif self.Server == "finished":
      print """It appears this sockets ring has finished.
To restart the ring, first delete the file server.txt"""
      return

    if self.IsServer:
      #WE ARE THE MAIN SERVER
      if self.Verbose: print "Starting sockets ring as main server"
      
      #initialize the arrays of clients
      self.Clients = []
      #start the server listening socket
      self.SSock = SSockClass()
      self.SSock.Listen(self.Port, self.MaxConnect)
      #update the sockets.txt file
      SocketsFile = os.path.join(self.Path, "sockets.txt")
      file(SocketsFile, "w").write("Initialized main server %s\n" % self.Server)
    
    else:
      #we are not the main server; setup as client
      if self.Verbose: print "Server is " + self.Server
      #create socket variable
      self.SSock = SSockClass()

    #flush std out
    if self.Verbose: sys.stdout.flush()


  def __GetServer(self):
    """Returns the name of the server and the port number,
or an empty string and zero if not found."""
    ServerFile = os.path.join(self.Path, "server.txt")
    if os.path.isfile(ServerFile):
      s = file(ServerFile).readline().strip().split()
      return s[0], int(s[1])
    else:
      return "", 0

  def __SetServer(self, ServerName, Port):
    ServerFile = os.path.join(self.Path, "server.txt")
    f = file(ServerFile, "w")
    f.write(ServerName + " " + str(Port) + "\n")
    f.close()

  def __DecideServer(self, PortMin = PortMinDflt, PortMax = PortMaxDflt):
    """Produces a file called server.txt in Path with
the name of a main server, which is automatically selected."""
    #first wait a random amt of time so two processors don't try to be server simultaneously
    w = random.random() * RandTime
    time.sleep(w)
    #get a port in case we need to be a server
    MyPort = GetPort(PortMin, PortMax)
    #check for another server
    ServerName, ServerPort = self.__GetServer()
    if ServerName == "":
      if MyPort == 0:
        #the port is bound so cannot be server
        self.Server = ""
        self.Port = 0
        self.IsServer = False
      else:
        #stake our claim as server and write name and port to file
        ServerName = socket.gethostname()
        self.Server = ServerName
        self.Port = MyPort
        self.__SetServer(ServerName, MyPort)
        self.IsServer = True
    else:
      #we are not the main server; setup as client
      self.Server = ServerName
      self.Port = ServerPort
      self.IsServer = False

  def Finished(self):
    """Indicates if the socketring has been finished."""
    return self.Server == "finished"


  #CLIENT MANIPULATION
      
  def ClientInd(self, s):
    "Returns the index of a particular socket."
    #see if s is a socket or ssockclass
    if type(s) is socket._socketobject:
      SockArray = [c.Sock for c in self.Clients]
      if s in SockArray:
        return SockArray.index(s)
      else:
        return None
    else:
      if s in self.Clients:
        return self.Clients.index(s)
      else:
        return None

  def GetClient(self, Sock):
    "Returns the SSockClass object for a given Sock."
    SockArray = [c for c in self.Clients if c.Sock is Sock]
    if len(SockArray) > 0:
      return SockArray[0]
    else:
      return None

  def __UpdateStats(self, Msg):
    """Updates the sockets.txt file."""
    SocketsFile = os.path.join(self.Path, "sockets.txt")
    if os.path.isfile(SocketsFile):
      #if the file is bigger than the max size, start over
      if os.path.getsize(SocketsFile) > MaxStatsSize:
        file(SocketsFile, "w").write("CONTINUED:\n" + Msg + "\n")
      else:
        #otherwise, append the message
        file(SocketsFile, "a").write(Msg + "\n")   

  def __RemoveClient(self, s, LogMsg = None):
    """Removes a socket from the client array."""
    #check if s is an index or the actual socket
    if type(s) is int:
      Ind = s
      SSock = self.Clients[Ind]
    else:
      Ind = self.ClientInd(s)
      SSock = s
    #check if there is an associated task; if so unlink it from the client
    if not SSock.Task is None:
      SSock.Task.Client = None
    if LogMsg is None:
      self.__UpdateStats("Deleting sock %s: %s (%d total)" % \
                         (SSock.Addr, SSock.Name, len(self.Clients)-1))
    else:
      self.__UpdateStats("Deleting sock %s: %s (%d total)  [%s]" % \
                         (SSock.Addr, SSock.Name, len(self.Clients)-1, LogMsg))      
    #close the socket
    SSock.Close()
    #delete the socket from the client array
    del self.Clients[Ind]

  def __AddClient(self, s):
    """Adds a socket client."""
    #check if s is a socket or a ssockclass
    if type(s) is socket._socketobject:
      SSock = SSockClass(s)
    else:
      SSock = s
    self.Clients.append(SSock)
    self.__UpdateStats("Adding sock %s: %s  (%d total)" % (SSock.Addr, SSock.Name, len(self.Clients)))

  def NextClient(self):
    """Returns the next non-busy client."""
    for c in self.Clients:
      if c.Connected and not c.Busy: return c
    return None

  def PruneClients(self):
    """Removes all closed clients."""
    self.Clients = [c for c in self.Clients if not c.Closed]

  def PruneHungClients(self, ClientTimeout = None):
    """Removes clients which have been running for too long."""
    if self.IsServer:
      if ClientTimeout is None: ClientTimeout = self.ClientTimeout
      i = 0
      while i < len(self.Clients):
        #get the client
        Client = self.Clients[i]
        if not Client.Busy:
          i += 1
          continue
        #see how long it's been running
        if not ClientTimeout is None \
          and Client.ElapsedSec() > ClientTimeout:
          #delete the client
          self.__RemoveClient(Client, 'exceeded client timeout')
        elif not Client.Task is None and not Client.Task.Timeout is None \
          and Client.ElapsedSec() > Client.Task.Timeout:
          #delete the client
          self.__RemoveClient(Client, 'exceeded task timeout')          
        else:
          i += 1

  def FreeClients(self):
    """Returns all non-busy clients."""
    return [x for x in self.Clients if not x.Busy]


  #TASK RUNNING ROUTINES

  def AddTask(self, TaskStr, Data = None, Timeout = None):
    """Adds a task to the task list and returns a task class."""
    self.LastID += 1
    NewTask = TaskClass(TaskStr, self.LastID, Data, Timeout)
    self.Tasks.append(NewTask)
    return NewTask

  def RemoveTask(self, Task):
    """Removes a task from the task list."""
    if Task in self.Tasks:
      Task.Client = None
      self.Tasks = [t for t in self.Tasks if not t is Task]

  def NextTask(self):
    """Returns the next task to run."""
    for t in self.Tasks:
      if t.Client is None and not t.Done: return t
    return None

  def GetRunTasks(self):
    """Returns all tasks currently running or to be run."""
    return [t for t in self.Tasks if not t.Done]

  def AllTasksDone(self):
    """Returns True if all tasks are complete."""
    return len(self.GetRunTasks()) == 0


  #SERVER ROUTINES

  def ListenForClients(self, TimeInSec = 10.):
    """Connects client sockets for use in tasks."""
    StartTime = time.time()
    while time.time() - StartTime < TimeInSec:
      #check for any socket requests
      SockArray = [self.SSock.Sock]
      (readers, writers, errors) = select.select(SockArray, [], SockArray, self.Interval)
      #check for requests
      for r in readers:
        #accept the request and add the socket
        (Sock, Addr) = r.accept()
        self.__AddClient(Sock)
      #check for sockets with errors
      for r in errors:
        #get the client
        Client = self.GetClient(r)
        #remove the client
        self.__RemoveClient(Client, 'error')
      #release unneeded clients
      self.ReleaseExcessClients()

  def RunTasks(self, TimeInSec = 10., ClientTimeout = None):
    """Runs as many tasks as possible from the task list."""
    StartTime = time.time()
    while (time.time() - StartTime < TimeInSec or TimeInSec == 0.):
      #try to send new tasks
      while True:
        #find the next non busy client and next non-running task
        Client = self.NextClient()
        Task = self.NextTask()
        #break if either is missing
        if Client is None or Task is None: break
        #assign the client and send the task
        if Client.SendExec(Task.TaskStr, Task.Data):
          #if the send went through, link the task and client
          Task.Client = Client
          Client.Task = Task
        else:
          #otherwise there was an error; remove this client
          self.__RemoveClient(Client, 'error')
      #check for any socket requests
      SockArray = [c.Sock for c in self.Clients] + [self.SSock.Sock]
      (readers, writers, errors) = select.select(SockArray, [], SockArray, self.Interval)
      #check for requests
      for r in readers:
        #see if another proc has joined
        if r is self.SSock.Sock:
          #accept the request and add the socket
          (Sock, Addr) = r.accept()
          self.__AddClient(Sock)
        else:
          #check the messages from existing procs; get the client and task
          Client = self.GetClient(r)
          Task = Client.Task
          #get the message
          Ret = Client.Recv()
          #did the processor finish?
          if not Ret is None:
            #mark the task as done
            Task.Done = True
            Task.TimeStop = time.time()
            #unlink the client
            Client.Task = None
            #check for a result
            if not Ret == MsgDone:
              Task.Result = Ret
          else:
            #something must have gone wrong;
            #remove the client
            self.__RemoveClient(Client, 'error')
      #check for sockets with errors
      for r in errors:
        #get the client
        Client = self.GetClient(r)
        #remove the client
        self.__RemoveClient(Client, 'error')
      #release unneeded clients
      self.ReleaseExcessClients()
      #release hung clients
      self.PruneHungClients(ClientTimeout)
      #if time is zero, we just want to do one loop
      if TimeInSec == 0.: break
      #if everything is done, we just want to do one loop
      if self.AllTasksDone(): break

  def ReleaseClients(self, Remove = False):
    """Stops all clients from waiting for the server for tasks."""
    if self.IsServer:
      i = 0
      while i < len(self.Clients):
        #get the client
        Client = self.Clients[i]
        #send a finish message
        if Client.Send(MsgRelease):
          if Remove:
            self.__RemoveClient(Client, Msg = 'released')
          else:
            i += 1
        else:
          #remove the socket
          self.__RemoveClient(Client, Msg = 'error')

  def ReleaseExcessClients(self, Remove = True):
    """Stops all clients from waiting for the server for tasks
    for those which exceed the number needed."""
    if self.IsServer:
      if self.MaxNeeded <= 0: return
      NExcess = len(self.Clients) - self.MaxNeeded
      for i in range(0, NExcess):
        #get the next non-busy client
        Client = self.NextClient()
        #make sure there are some
        if Client is None: return
        #send a finish message
        if not Client.Send(MsgRelease):
          #error; remove the socket
          self.__RemoveClient(Client, 'error')
        elif Remove:
          self.__RemoveClient(Client, 'released as excess')


  def Efficiency(self):
    """Gives the overall efficiency in terms of fraction client idle time."""
    TotWork = sum([c.SecWork for c in self.Clients])
    TotIdle = sum([c.SecIdle for c in self.Clients])
    return TotWork / (TotWork + TotIdle + 0.001)


  #BOTH SERVER & CLIENT

  def Finish(self):
    """Closes down all clients in the socket ring."""
    if self.IsServer:
      while len(self.Clients) > 0:
        #get the last client
        Client = self.Clients[-1]
        #send a finish message
        Client.Send(MsgFinish)
        #remove the socket
        self.__RemoveClient(Client, 'finished')
      #close down the listening socket
      self.SSock.Close()
      self.SSock = None
      #update the server.txt file
      self.__SetServer("finished", 0)
      self.Server = "finished"
      self.Port = 0
      if self.Verbose: print "Finished running as socket main server."
    else:
      self.SSock.Close()
      self.Server = "finished"
      self.Port = 0
      

  #CLIENT ROUTINES

  def ConnectServer(self, Timeout = None):
    """Connects the client to the server."""
    #try to connect to the server
    Attempts = 0
    if Timeout is None: Timeout = self.MaxIdle
    StartTime = time.time()
    while time.time() - StartTime < Timeout:
      Attempts += 1
      if Attempts < 10 and self.Verbose:
        print "Connecting to server %s on port %d" % (self.Server, self.Port)
      #connect to the socket
      if self.SSock.Connect(self.Server, self.Port):
        print "Connected!"
        return True
      else:
        if Attempts < 10 and self.Verbose:
          print "Socket error connecting to server %s on port %d" % (self.Server, self.Port)
        #sleep a short interval before retrying
        time.sleep(self.Interval)
    return False

    
  def RunClient(self):
    """Runs the client and sends it tasks."""
    #check if the simulation has finished
    if self.Server == "finished":
      return True
    #set the time of the last non-idle event
    StartTime = time.time()
    Cont = True
    if self.Verbose: print "Running as socket client"
    while Cont and time.time() - StartTime < self.MaxIdle:
      if self.SSock.Connected:
        #get the next command
        Msg = self.SSock.Recv()
        StartTime = time.time()
        #see if this is a finish message
        if Msg == MsgFinish:
          #stop trying to connect
          Cont = False
          #close down the socket
          self.SSock.Close()
          if self.Verbose:
            print "\nEnded connection with server %s on port %d" % (self.Server, self.Port)
            print "Finished running as socket client"
          return True
        elif Msg == MsgRelease:
          #stop trying to connect
          Cont = False
          if self.Verbose:
            print "\nStopped running tasks from server %s on port %d" % (self.Server, self.Port)
          return True
        elif Msg is None:
          #there was an error
          self.SSock.Close()
          if self.Verbose:
            print "\nDisconnected from server %s on port %d" % (self.Server, self.Port)
        else:
          #it was some other command, print a dot to indicate we did it
          if self.Verbose and TaskSignal:
            sys.stdout.write(".")
            NWrite += 1
            if NWrite >= 70:
              sys.stdout.write("\n")
              NWrite = 0
      else:
        #if the sock is closed, make a new connection
        if self.SSock.Closed:
          self.SSock.GetNewSock()
        #try to connect to the server
        if self.ConnectServer():
          #zero the number of outputs
          NWrite = 0
          #update the last event time
          StartTime = time.time()
    if self.Verbose:
      print "Finished running as socket client."
    #we didn't receive a finish message
    return False


