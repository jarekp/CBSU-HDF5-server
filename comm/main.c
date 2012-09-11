/* 
 * File:   main.c comm 
 * Author: jarekp
 *
 * Created on February 6, 2012, 3:03 PM
 */

#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <time.h>
#include <string.h>
#include <crypt.h>
#include <fcntl.h>
#include <errno.h>
#include <unistd.h>
#include <signal.h>

#define MAXLINE 5000000

#define INSTALLDIR "/usr/local/hdf5"

char prtbuf[MAXLINE];
char prtbuf1[500];
char lockfile[1000];
char hdf5dir[256], hdf5conf[256], tmpdir[256], hdf5root[512];
char *remote_host;

char **cqueue;
int nc;
int pid;

FILE *flog;
FILE *outs, *ins;

void chomp(char *inp)
{
    if (inp[0] == '\0')return;
    if (inp[strlen(inp) - 1] == '\n')inp[strlen(inp) - 1] = '\0';
}

void logmsg(char *msg)
{
    char buf[100];

    time_t ttt = time(NULL);
    sprintf(buf, "%s", asctime(localtime(&ttt)));
    chomp(buf);
    fprintf(flog, "%d %s [%s]: %s", pid, remote_host, buf, msg);
    fflush(flog);
}

void error(char *msgin, char *msgout)
{
    logmsg(msgin);
    fprintf(stdout, "%s\n", msgout);
    fflush(stdout);
    fclose(flog);
    exit(1);
}

void error1(char *msgin, char *msgout, int line)
{
    char buf[1000];

    sprintf(buf, "%d %s\n", line, msgin);
    error(buf, msgout);
}

void printout(char *msg)
{
    logmsg(msg);
    fprintf(stdout, msg);
    fflush(stdout);
}

int getlock(int sleep)
{
    struct stat stt;
    FILE *lll;
    time_t ttt;
    clock_t c_start, c_end;

    int i = stat(lockfile, &stt);
    if (i != 0)
    {
        //no lock file - getting lock
        if ((lll = fopen(lockfile, "w")) == NULL)
        {
            return -2;
        }
        ttt = time(NULL);
        fprintf(lll, "%s\n", asctime(localtime(&ttt)));
        fclose(lll);
        return 1;
    }
    else
    {
        //lock file is there
        //sleep for sleep/1000 seconds
        c_start = clock();
        while (1)
        {
            c_end = clock();
            double etime = (double) (c_end - c_start) / CLOCKS_PER_SEC;
            if (etime * 1000 > (double) sleep)break;
        }
        return -1;
    }

    return 1;
}

int main(int argc, char** argv)
{

    char command[50];
    char buf[1000];
    FILE *in;
    time_t ttt;

    strcpy(hdf5conf, INSTALLDIR);

    pid = getpid();
    sprintf(prtbuf, "%s/cbsuhdf5comm.log", hdf5conf);
    if ((flog = fopen(prtbuf, "a")) == NULL)
    {
        fprintf(stderr, "ERROR: Cannot open log file %s\n", prtbuf);
        exit(1);
    }

    remote_host = getenv("REMOTE_HOST");
    sprintf(prtbuf, "%s/config.txt", hdf5conf);
    //read config file first
    if ((in = fopen(prtbuf, "r")) != NULL)
    {
        if (!fgets(prtbuf, 4000, in))error1("EOF", "Service configuration error 1", __LINE__);
        if (!fgets(prtbuf, 100, in))error1("EOF", "Service configuration error 2", __LINE__);
        if (!fgets(prtbuf, 4000, in))error1("EOF", "Service configuration error 3", __LINE__);
        if (!fgets(hdf5dir, 500, in))error1("EOF", "Service configuration error 4", __LINE__);
        chomp(hdf5dir);
        ttt = time(NULL);
        logmsg("starting\n");
        if (!fgets(prtbuf, 4000, in))error1("EOF", "Service configuration error 5", __LINE__);
        if (!fgets(hdf5root, 500, in))error1("EOF", "Service configuration error 6", __LINE__);
        chomp(hdf5root);
        sprintf(prtbuf, "%s%s", hdf5dir, hdf5root);
        strcpy(hdf5root, prtbuf);
        if (!fgets(prtbuf, 4000, in))error1("EOF", "Service configuration error 7", __LINE__);
        if (!fgets(prtbuf, 500, in))error1("EOF", "Service configuration error 8", __LINE__);
        if (!fgets(prtbuf, 4000, in))error1("EOF", "Service configuration error 9", __LINE__);
        if (!fgets(command, 50, in))error1("EOF", "Service configuration error 10", __LINE__);
        chomp(command);
        fclose(in);
        sprintf(lockfile, "%s/hdf5.lock", hdf5dir);

        //check if the server is alive
        sprintf(prtbuf, "%s/cbsuhdf5.pid", hdf5conf);
        int restart = 0;
        char pid[40];
        if ((in = fopen(prtbuf, "r")) != NULL)
        {
            if (fgets(pid, 39, in))
            {
                chomp(pid);
                sprintf(prtbuf, "/proc/%s/cmdline", pid);
                fclose(in);
                if ((in = fopen(prtbuf, "r")) != NULL)
                {
                    restart = 0;
                    fclose(in);
                }
                else
                {
                    restart = 1;
                }
            }
            else
            {
                fclose(in);
                restart = 1;
            }
        }
        else
        {
            restart = 1;
        }
        if (restart == 1)
        {
            logmsg("restarting hdf5 demon!\n");
            sprintf(prtbuf, "rm -f %s/cbsuhdf5.pid %s/cbsuhdf5.pid1 1> /dev/null 2> /dev/null", hdf5conf, hdf5conf);
            int jk = system(prtbuf);
            sprintf(prtbuf, "%s/cbsuhdf5d start 1> /dev/null 2> /dev/null", hdf5conf);
            jk = system(prtbuf);
            logmsg(" .. done\n");
        }

        //set communication mode
        if (command[0] == '_')
        {
            sprintf(prtbuf, "%s%s_out", hdf5dir, &command[1]);
            sprintf(prtbuf1, "Opening %s pipe for reading\n", prtbuf);
            logmsg(prtbuf1);
            if ((outs = fopen(prtbuf, "r")) == NULL)
            {
                sprintf(buf, "Cannot open output pipe descriptor %s\n", prtbuf);
                error(buf, "Service error 1");
            }
            sprintf(prtbuf, "%s%s_in", hdf5dir, &command[1]);
            sprintf(prtbuf1, "Opening %s pipe for writing\n", prtbuf);
            logmsg(prtbuf1);
            if ((ins = fopen(prtbuf, "w")) == NULL)
            {
                sprintf(buf, "Cannot open input pipe descriptor %s\n", prtbuf);
                error(buf, "Service error 2");
            }
        }
        else
        {
            //invalid option
            //nothing to do anyway
            error("ERROR: No server communicator defined\n", "Service error 3");
        }
    }
    else
    {
        error("ERROR: Cannot open config.txt\n", "Service input error 1");
    }
    //first get the input and store it in memory
    nc = 0;
    char cmd[50];
    while (1)
    {
        int i = 0;
        while (1)
        {
            char c = fgetc(stdin);
            if (c != 13)
            {
                prtbuf[i] = c;
                i++;
            }
            if (c == 10)
            {
                prtbuf[i] = '\0';
                break;
            }
        }
        //fprintf(stdout, "==%s==", prtbuf);
        if (prtbuf[0] == '\n')
        {
            if (strcmp(cmd, "QUERY\n") == 0 && (nc == 5 || nc == 9))
            {
                //ignore
            }
            else if (strcmp(cmd, "TABLE\n") == 0 && nc == 5)
            {
                //ignore
            }
            else
            {
                break;
            }
        }
        if (nc == 0)
        {
            cqueue = malloc(sizeof (char *));
            if (cqueue == NULL)error("Cannot allocate memory (1)\n", "Service memory error 1");
        }
        else
        {
            cqueue = realloc(cqueue, (nc + 1) * sizeof (char *));
            if (cqueue == NULL)error("Cannot allocate memory (2)\n", "Service memory error 2");
        }
        cqueue[nc] = malloc(strlen(prtbuf) + 1);
        if (cqueue[nc] == NULL)error("Cannot allocate memory (3)\n", "Service memory error 3");
        strcpy(cqueue[nc], prtbuf);
        if (nc == 0)
        {
            strcpy(cmd, prtbuf);
        }
        nc++;
    }
    //get a lock on the server and push in the input
    sprintf(buf, "%d lines received\n", nc);
    logmsg(buf);
    int ii = 0;
    int iret = 0;
    while ((iret = getlock(100)) < 0)
    {
        if (iret == -2)
        {
            error("ERROR: cannot create a lock\n", "Service lock error 1");
            exit(1);
        }
        printout("server busy, waiting\n");
        ii++;
        if (ii == 100)error("ERROR: connection timeout\n", "server timeout error");
    }
    sprintf(buf, "lock established (%d)\n", ii);
    logmsg(buf);
    logmsg(buf);
    //push the input to the server
    for (ii = 0; ii < nc; ii++)
    {
        fprintf(ins, "%s", cqueue[ii]);
        fflush(ins);
        free(cqueue[ii]);
    }
    fprintf(ins, "\n");
    fflush(ins);
    logmsg("data fed to server\n");
    //transmit output to the client
    while (1)
    {
        fgets(prtbuf, MAXLINE - 1, outs);
        fprintf(stdout, "%s", prtbuf);
        fflush(stdout);
        if (strcmp(prtbuf, "Command COMPLETED\n") == 0)break;
    }
    logmsg("Command COMPLETED\n");
    //remove the lock
    iret = unlink(lockfile);
    if (iret == -1)
    {
        logmsg("Cannot remove lock\n");
        sprintf(prtbuf, "rm -f %s", lockfile);
        system(prtbuf);
    }
    //free memory
    logmsg("done\n");
    free(cqueue);
}

