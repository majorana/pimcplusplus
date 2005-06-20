#ifndef RUN_INFO_H
#define RUN_INFO_H

#include "Common/IO/InputOutput.h"
#include <unistd.h>
#include <stdio.h>
#include <sys/types.h>
#include <pwd.h>
#include <time.h>
class RunInfoClass
{
public:
  string ProgramName;
  string Version;
  string CommonVersion;
  string UserName;
  string RunTime;
  string BuildDate;
  string BuildTime;
  string HostName;
  inline void Write(IOSectionClass &outSec)
  {
    outSec.WriteVar("ProgramName", ProgramName);
    outSec.WriteVar("Version", Version);
    outSec.WriteVar("CommonVersion", CommonVersion);
    outSec.WriteVar("UserName", UserName);
    outSec.WriteVar("RunTime", RunTime);
    outSec.WriteVar("BuildDate", BuildDate);
    outSec.WriteVar("BuildTime", BuildTime);
    outSec.WriteVar("HostName", HostName);
  }
  inline void Read(IOSectionClass &inSec)
  {
    assert (inSec.ReadVar("ProgramName", ProgramName));
    assert (inSec.ReadVar("Version", Version));
    if (!inSec.ReadVar ("CommonVersion", CommonVersion))
      CommonVersion = "unknown";
    assert (inSec.ReadVar("UserName", UserName));
    assert (inSec.ReadVar("RunTime", RunTime));
    assert (inSec.ReadVar("BuildDate", BuildDate));
    assert (inSec.ReadVar("BuildTime", BuildTime));
    assert (inSec.ReadVar("HostName", HostName));
  }
  inline RunInfoClass()
  {
    struct passwd* pwInfo = getpwuid(getuid());
    UserName = pwInfo->pw_name;
    BuildDate = __DATE__;
    BuildTime = __TIME__;
    char hostname[300];
    gethostname(hostname, 300);
    HostName = hostname;
    time_t seconds = time(NULL);
    RunTime = ctime(&seconds);
    Version = VERSION;
    CommonVersion = COMMONVERSION;
  }
};


#endif
