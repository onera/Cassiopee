/*    
    Copyright 2013-2019 Onera.

    This file is part of Cassiopee.

    Cassiopee is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Cassiopee is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Cassiopee.  If not, see <http://www.gnu.org/licenses/>.
*/

// cplot_server support

#include "GenIO.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <sys/types.h>
#if defined(_WIN32) ||defined(_WIN64)
# include <winsock2.h>
# include <ws2tcpip.h>
#else
# include <sys/socket.h>
# include <netinet/in.h>
# include <netdb.h>
# include <unistd.h>
# include <sys/utsname.h>
# include <sys/time.h>
#endif
#include <signal.h>
#include <fcntl.h>
#include <errno.h>
#include <time.h>

//=============================================================================
/* cplotSendDouble
   Send a double to cplot_server */
//=============================================================================
void K_IO::GenIO::cplotSendDouble(int* mySocket, 
                                  double value)
{
  int ndouble = sizeof(double);
  
  if (write(*mySocket, &value, ndouble) < 0)
  {
    perror("write");
    exit(1);
  }
}

//=============================================================================
/* cplotSendInteger
   Send an integer to cplot_server */
//=============================================================================
void K_IO::GenIO::cplotSendInteger(int* mySocket, 
                                   int value)
{
  int nint = sizeof(int);

  if (write(*mySocket, &value, nint) < 0)
  {
    perror("write");
    exit(1);
  }
}

//=============================================================================
/* cplotSendHeader
   Send the header of the grid to cplot_server */
//=============================================================================
void K_IO::GenIO::cplotSendHeader(int* mySocket, int nvar, int nfield, 
                                  int zonenumber, int ni, int nj, int nk)
{
  int nint = sizeof(int);
  
  if (write(*mySocket, &nvar, nint) < 0)
  {
    perror("write");
    exit(1);
  }

  if (write(*mySocket, &nfield, nint) < 0)
  {
    perror("write");
    exit(1);
  }

  if (write(*mySocket, &zonenumber, nint) < 0)
  {
    perror("write");
    exit(1);
  }

  if (write(*mySocket, &ni, nint) < 0)
  {
    perror("write");
    exit(1);
  }

  if (write(*mySocket, &nj, nint) < 0)
  {
    perror("write");
    exit(1);
  }

  if (write(*mySocket, &nk, nint) < 0)
  {
    perror("write");
    exit(1);
  }
}

//=============================================================================
/* cplotSend
   Send array of doubles of size n to server */
//=============================================================================
void K_IO::GenIO::cplotSend(int* mySocket, int size, double* array)
{
  int ndouble = sizeof(double);
  
  if (write(*mySocket, array, size*ndouble) < 0)
  {
    perror("write");
    exit(1);
  }
}

//=============================================================================
/* Cree la socket pour le client */
//=============================================================================
E_Int K_IO::GenIO::cplotClient(char* machine, char* service, int *mySocket)
{
  struct hostent *hostp;
  struct servent *servp;
  struct sockaddr_in server;
  unsigned int sock;
  fd_set mask;

  if ((sock = socket(AF_INET, SOCK_STREAM, IPPROTO_TCP)) < 0)
  {
    perror("socket");
    return 1;
  }
  
  if (isdigit(service[0]))
  {
    static struct servent s;
    servp = &s;
    s.s_port = htons((u_short)atoi(service));
  }
  else if ((servp = getservbyname(service, "tcp")) == 0)
  {
    printf("connecting cplot: service=%s on machine=%s\n", service, machine);
    fprintf(stderr, "%s: unknown service.\n", service);
    return 1;
  }
  if ((hostp = gethostbyname(machine)) == 0)
  {
    printf("connecting cplot: service=%s on machine=%s\n", service, machine);
    fprintf(stderr, "%s: unknown host.\n", machine);
    return 1;
  }
  memset((void*) &server, 0, sizeof server);

  server.sin_family = AF_INET;
  memcpy((void*) &server.sin_addr, hostp->h_addr, hostp->h_length);
  server.sin_port = servp->s_port;
  if (connect(sock, (struct sockaddr *)&server, sizeof server) < 0)
  {
    (void) close(sock);
    perror("connect");
    return 1;
  }
  FD_ZERO(&mask);
  FD_SET(sock, &mask);
  FD_SET(fileno(stdin), &mask);

  *mySocket = sock;
  printf("connecting cplot: service=%s on machine=%s\n", service, machine);
  return 0;
}
