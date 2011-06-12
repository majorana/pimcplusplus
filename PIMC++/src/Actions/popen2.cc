#include "popen2.h"
#include <cstdlib>

pid_t popen2(const char *shell_cmd, int *p_fd_in, int *p_fd_out)
{
  //CREATING TWO PIPES:
  int fds_processInput[2];  //pipe for process input
  int fds_processOutput[2]; //pipe for process output
  
  if(pipe(fds_processInput) != 0) //create process input pipe
    {
      std::cerr << "pipe (process input) failed\n";
      exit(1);
    }
  
  if(pipe(fds_processOutput) != 0) //create process output pipe
    {
      std::cerr << "pipe (process output) failed\n";
      exit(1);
    }
  
  //FORKING A CHILD PROCESS:
  pid_t pid;
  if((pid = fork()) < 0)
    {
      std::cerr << "fork failed\n";
      exit(2);
    }
   
  //CONNECT THE CORRECT PIPE ENDS IN THE CHILD:
  if(pid == 0)  //child process
    {
      //for process input pipe:
      close(fds_processInput[1]);   //close output
      dup2(fds_processInput[0], 0); //close fd 0, fd 0 = fds_processInput[0]  
      //for process output pipe:
      close(fds_processOutput[0]);   //close input
      dup2(fds_processOutput[1], 1); //close fd 1, fd 1 = fds_processOutput[1]
      
      
      execl("/bin/sh", "sh", "-c", shell_cmd, 0 ); 
      std::cerr << "failed to run shell_cmd\n";
    }
  else  //parent process
    {
      //for process input pipe:
      close(fds_processInput[0]);   //close input
       
      //for process output pipe:
      close(fds_processOutput[1]);   //close output

      if(p_fd_in == 0)
	close(fds_processInput[1]);
      else
	*p_fd_in = fds_processInput[1];
       
      if(p_fd_out == 0)
	close(fds_processOutput[0]);
      else
	*p_fd_out = fds_processOutput[0];

    }
  return pid; 
}

std::string myRead(int fd_read){
  std::string tag = "";
  char c;
  read(fd_read, &c, 1);
  while(c != '\n'){
    tag += c;
    read(fd_read, &c, 1);
  }
  return tag;
}

int myWrite(int fd_write, std::string send){
  send += '\n';
  //int size = send.size() + 1;
  int size = send.size();
  return(write(fd_write, send.c_str(), size));
}




//int main(int argc, char** argv)
//{
//  MPI_Init(&argc, &argv);
//  
//  int rank, size;
//  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//  MPI_Comm_size(MPI_COMM_WORLD, &size);
//  
//  cout << "Hello World! I am " << rank << " of " << size << endl;
//  if(rank == 0){
//    int toread;
//    int towrite;
//    std::string tag;
//    cerr << rank << " Attempting popen...";
//    popen2("mpirun -np 1 ./comp",&towrite,&toread);
//    cerr << " Successful" << endl;
//
//    //write(towrite,"INIT\n",5);
//    //read(toread,line,16);
//    cerr << rank << " attempting read: ";
//    tag = myRead(toread);
//    cerr << "I picked up " << tag << endl;
// 
//    ///reads the input or prints no data yet if nothign there...replace this with while (I don't have the full input I want yet (including the '\n')
//
//    int send = write(towrite,"INIT\n",5);
//    while (1){
//      //read(toread,line,16);
//      tag = myRead(toread);
//      cerr << "I picked up " << tag << endl;
//      if(tag == "DATA"){
//        //read(toread,line,16);
//        tag = myRead(toread);
//        cerr << "I picked up data " << tag << endl;
//        int send = write(towrite,"INIT\n",5);
//        cerr << "Run sent " << send << endl;
//      //if (read(toread,&c,1)<0)
//      //  cerr<<"No data yet"<<endl;
//      //else {
//      //  cerr<<c;
//      //}
//      }
//    }
//  }
//  MPI_Finalize();
//}
