#ifndef COMMUNICATOR_CLASS_H
#define COMMUNICATOR_CLASS_H

class CommunicatorClass
{
 public:
  int NumProcs() const;
  int MyProc() const;

  void SendReceive (int SendProc, const Array<dVec,1> &SendBuff,
		    int RecvProc,       Array<dVec,1> &RecvBuff);

  void SendReceive (int SendProc, const Array<double,1> &SendBuff,
		    int RecvProc,       Array<dVec,1> &RecvBuff);

  void SendReceive (int SendProc, const Array<int,1> &SendBuff,
		    int RecvProc,       Array<dVec,1> &RecvBuff);
};



#endif
