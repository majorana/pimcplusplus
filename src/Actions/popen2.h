#ifndef POPEN_2_UTILS_H
#define POPEN_2_UTILS_H

#include <stdio.h>
#include <sys/types.h>
#include <unistd.h>
#include <iostream.h>
#include <istream.h>
#include <string>

pid_t popen2(const char *shell_cmd, int *p_fd_in, int *p_fd_out);

std::string myRead(int fd_read);

int myWrite(int fd_write, std::string send);
  
#endif
