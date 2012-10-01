/* 
 * File:   main.c
 * Author: jarekp
 *
 * Created on December 20, 2011, 3:29 PM
 */

#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <time.h>
#include <string.h>
#include <crypt.h>
#include <fcntl.h>
#include <errno.h>
#include <unistd.h>
#include <signal.h>
#include "hdf5.h"

/*
 * 
 */

#define ITERLIMIT  100

#define DATAARRAY "data_array"
#define MAXSTRLEN  256

#define INSTALLDIR "/usr/local/hdf5"

typedef struct
{
    char username[50];
    char password[50];
} userdata_type;

typedef struct
{
    char *string;
    long unsigned int index;
} arrsort;

typedef struct
{
    char flags; // bit 1 strand, bit 2-3 allele, bit 4-8 undef
    char nums[8];
} snpmetadata;

typedef struct
{
    char name[256];
    int num;
    int mnt;
    int order; //1 - positions fast dim, 2 - taxa fast dim, 3 - both present
    int enc;
    char encoding[50];
} projects_type;

typedef struct
{
    char name[256];
    int mnt;
} files_type;

typedef struct
{
    char name[11];
    int mnt;
    int prj;
    int num;
    hsize_t dim0; //positions dim
    hsize_t dim1; //taxa dim
} chrinfo_type;

typedef struct
{
    char cell[3];
} threebyte_type;

char pass[100];
char tmpdir[500];
char hdf5root[500];
char hdf5dir[500];
char hdf5conf[500];
char prtbuf[5000];
char errbuf[5000];
char *outsbuffer;
int delayed;
FILE *flog;
FILE *outs;
FILE *ins;
int i_outs, i_ins;
hid_t file;
projects_type *projects;
int nprojects;
files_type *files;
int nfiles;
hid_t *filesh;
int nchrinfo;
chrinfo_type *chrinfo;
int comm_pid;
struct timeval c_start, c_start0, c_end;
int nusers;
userdata_type *users;
char **accesstable; //same order as users: '\npname1\npname2\npname3\n'


int query(), mount(), umount(), pinfo(), plist(), flist(), finfo(), cache(), quit(), sort(), table(), useradd(), userdel(), userpass(), useracc(), aflist(), ulist(), login();
int read_access_table(), write_access_table(), umount1(), delproject(), printtable();
void freestrtable();

void shutdown();

void chomp(char *inp)
{
    if (inp[0] == '\0')return;
    if (inp[strlen(inp) - 1] == '\n')inp[strlen(inp) - 1] = '\0';
}

void printout(char *txt)
{
    fprintf(flog, "%s", txt);
    fflush(flog);
    if (delayed == 0)
    {
        fprintf(outs, "%s", txt);
        fflush(outs);
    }
    else
    {
        if (outsbuffer == NULL)
        {
            outsbuffer = malloc(strlen(txt) + 1);
            sprintf(outsbuffer, "%s", txt);
        }
        else
        {
            unsigned int pos = strlen(outsbuffer);
            outsbuffer = realloc(outsbuffer, pos + strlen(txt) + 1);
            strcat(outsbuffer, txt);
        }
    }
}

void printoutbuf()
{
    if (outsbuffer != NULL)
    {
        fprintf(outs, "%s", outsbuffer);
        fflush(outs);
        free(outsbuffer);
        outsbuffer = NULL;
    }
}

void error(char *msg, int line)
{
    fprintf(flog, "%d: %s\n", line, msg);
    fprintf(stderr, "%d: %s\n", line, msg);
    shutdown();
    fclose(flog);
    exit(1);
}

int readline(char *buf, int size, int line, int chmp)
{
    time_t tt;

    if (!fgets(buf, size, ins))
    {
        fprintf(flog, "%d: EOF\n", line);
        fprintf(stderr, "%d: EOF\n", line);
        fclose(flog);
        shutdown();
        exit(1);
    }
    if (buf[0] == '\n' && chmp != -1)
    {
        fprintf(outs, "end of input - executing\n");
        fprintf(flog, "end of input - executing\n");
        tt = time(NULL);
        fprintf(outs, "START %s", asctime(localtime(&tt)));
        fprintf(flog, "START %s", asctime(localtime(&tt)));
        fflush(flog);
        fflush(outs);
        gettimeofday(&c_start, NULL);
        return -1;
    }
    if (chmp == 1)
    {
        if (buf[strlen(buf) - 1] == '\n')buf[strlen(buf) - 1] = '\0';
    }
    gettimeofday(&c_start, NULL);
    fprintf(flog, "%s", buf);
    if (chmp == 1)fprintf(flog, "\n");
    fflush(flog);
    return 1;
}

void gotoend()
{
    while (fgets(prtbuf, 4999, ins))
    {
        if (prtbuf[0] == '\n')break;
    }
    gettimeofday(&c_start, NULL);
}

int find_project(char *name)
{
    int i;

    if (projects[0].mnt == -1 || nprojects == 0)return -1;
    for (i = 0; i < nprojects; i++)
    {
        //fprintf(flog, "%d\n", i);
        //fflush(flog);
        if (strcmp(projects[i].name, name) == 0)return i;
    }
    return -1;
}

int get_enc_size(char *buf)
{
    char onebyte[] = "|IUPAC|2hex_phased|2hex_unphased|AB_phased|AB_unphased|";
    char threebyte[] = "|2hex+_phased|2hex+_unphased|";

    if (strstr(onebyte, buf) == NULL)
    {
        if (strstr(threebyte, buf) == NULL)
        {
            return -1;
        }
        else
        {
            return 3;
        }
    }
    else
    {
        return 1;
    }
}

int find_user(char *name)
{
    int i;

    if (users[0].username[0] == '\0' || nusers == 0)return -1;
    for (i = 0; i < nusers; i++)
    {
        if (strcmp(users[i].username, name) == 0)return i;
    }
    return -1;
}

int find_file(char *name)
{
    int i;

    if (files[0].mnt == -1 || nfiles == 0)return -1;
    for (i = 0; i < nfiles; i++)
    {
        if (strcmp(files[i].name, name) == 0)return i;
    }
    return -1;
}

int find_chr(char *name, int prj)
{
    int i;

    if (nchrinfo == 0)return -1;
    for (i = 0; i < nchrinfo; i++)
    {
        if (strcmp(chrinfo[i].name, name) == 0 && chrinfo[i].prj == prj)return i;
    }
    return -1;
}

int read_str_attr(hid_t project, char *attr_name, char *attr_value)
{
    hid_t attr, atype;
    herr_t ret, ret1;
    H5A_info_t ainfo;

    if (!H5Aexists(project, attr_name))return -1;
    attr = H5Aopen_name(project, attr_name);
    H5Aget_info(attr, &ainfo);
    atype = H5Tcopy(H5T_C_S1);
    H5Tset_size(atype, ainfo.data_size);
    ret = H5Aread(attr, atype, attr_value);
    attr_value[ainfo.data_size] = '\0';
    ret1 = H5Aclose(attr);
    ret1 = H5Tclose(atype);
    if (ret < 0)return -1;
    return 1;
}

int read_int_attr(hid_t project, char *attr_name, int *attr_value)
{
    hid_t attr;
    herr_t ret;

    if (!H5Aexists(project, attr_name))return -1;
    attr = H5Aopen_name(project, attr_name);
    ret = H5Aread(attr, H5T_NATIVE_INT, attr_value);
    if (ret < 0)return -1;
    ret = H5Aclose(attr);
    return 1;
}

int read_chrinfo(hid_t project, int mnt, int prj, int test, int prt, int order)
{
    int i, n;
    hsize_t dimspf[2], dimstf[2], dimss[2], taxadim, markersdim, positionsdim;
    char buf[1000], orderstr[20];
    hid_t chromosome, dataset, filespace;
    herr_t ret;


    dataset = H5Dopen(project, "taxa", H5P_DEFAULT);
    if (dataset < 0)
    {
        sprintf(buf, "Cannot open dataset 'taxa'\n");
        printout(buf);
        return -1;
    }
    filespace = H5Dget_space(dataset);
    ret = H5Sget_simple_extent_dims(filespace, dimss, NULL);
    taxadim = dimss[0];
    H5Dclose(dataset);
    i = 0;
    n = 0;
    strcpy(orderstr, "positions-fast");
    if (order == 2)strcpy(orderstr, "taxa-fast");
    while (1)
    {
        n++;
        sprintf(buf, "chr%d", n);
        if (!H5Lexists(project, buf, H5P_DEFAULT))
        {
            break;
        }
        chromosome = H5Gopen(project, buf, H5P_DEFAULT);
        if (order == 1 || order == 3)
        {
            sprintf(buf, "%s_pf", DATAARRAY);
            dataset = H5Dopen(chromosome, buf, H5P_DEFAULT);
            if (dataset < 0)
            {
                sprintf(prtbuf, "Cannot open dataset '%s' in chr %d\n", buf, n);
                printout(prtbuf);
                return -1;
            }
            filespace = H5Dget_space(dataset);
            ret = H5Sget_simple_extent_dims(filespace, dimspf, NULL);
            H5Dclose(dataset);
            H5Sclose(filespace);
            if (dimspf[0] != taxadim)
            {
                sprintf(buf, "dimensions of '%s_pf' in project %d chr %d and taxa array don't match in declared order %s (%lli %lli %lli)\n", DATAARRAY, prj, n, orderstr, dimspf[0], dimspf[1], taxadim);
                printout(buf);
                return -1;
            }
        }
        if (order == 2 || order == 3)
        {
            sprintf(buf, "%s_tf", DATAARRAY);
            dataset = H5Dopen(chromosome, buf, H5P_DEFAULT);
            if (dataset < 0)
            {
                sprintf(prtbuf, "Cannot open dataset '%s' in chr %d\n", buf, n);
                printout(prtbuf);
                return -1;
            }
            filespace = H5Dget_space(dataset);
            ret = H5Sget_simple_extent_dims(filespace, dimstf, NULL);
            H5Dclose(dataset);
            H5Sclose(filespace);
            if (dimstf[1] != taxadim)
            {
                sprintf(buf, "dimensions of '%s_tf' in project %d chr %d and taxa array don't match in declared order %s (%lli %lli %lli)\n", DATAARRAY, prj, n, orderstr, dimstf[0], dimstf[1], taxadim);
                printout(buf);
                return -1;
            }
        }
        dataset = H5Dopen(chromosome, "positions", H5P_DEFAULT);
        if (dataset < 0)
        {
            sprintf(buf, "Cannot open dataset 'positions' in chr %d\n", n);
            printout(buf);
            return -1;
        }
        filespace = H5Dget_space(dataset);
        ret = H5Sget_simple_extent_dims(filespace, dimss, NULL);
        positionsdim = dimss[0];
        H5Dclose(dataset);
        H5Sclose(filespace);
        if (order == 1 || order == 3)
        {
            if (dimspf[1] != positionsdim)
            {
                sprintf(buf, "dimensions of '%s_pf' in prj %d chr %d and positions array don't match (%lli %lli)\n", DATAARRAY, prj, n, dimspf[1], positionsdim);
                printout(buf);
                return -1;
            }
        }
        if (order == 2 || order == 3)
        {
            if (dimstf[0] != positionsdim)
            {
                sprintf(buf, "dimensions of '%s_tf' in prj %d chr %d and positions array don't match (%lli %lli)\n", DATAARRAY, prj, n, dimstf[0], positionsdim);
                printout(buf);
                return -1;
            }
        }
        if (!H5Lexists(chromosome, "markers", H5P_DEFAULT))
        {
            sprintf(buf, "Cannot open dataset 'markers' in chr %d - markers not present (OK, they are optional)\n", n);
            printout(buf);
        }
        else
        {
            if (!H5Lexists(chromosome, "markers_sv", H5P_DEFAULT) || !H5Lexists(chromosome, "markers_si", H5P_DEFAULT))
            {
                if (prt == 1)
                {
                    sprintf(buf, "WARNING: 'markers' table of chr %d is not indexed - text searches are disabled\n", n);
                    printout(buf);
                }
            }
        }
        dataset = H5Dopen(chromosome, "alleles", H5P_DEFAULT);
        if (dataset < 0)
        {
            sprintf(buf, "Cannot open dataset 'alleles' in chr %d\n", n);
            printout(buf);
            return -1;
        }
        filespace = H5Dget_space(dataset);
        ret = H5Sget_simple_extent_dims(filespace, dimss, NULL);
        markersdim = dimss[0];
        H5Dclose(dataset);
        H5Sclose(filespace);
        if (order == 1 || order == 3)
        {
            if (dimspf[1] != markersdim)
            {
                sprintf(buf, "dimensions of '%s_pf' in prj %d chr %d and alleles array don't match (%lli %lli)\n", DATAARRAY, prj, n, dimspf[1], markersdim);
                printout(buf);
                return -1;
            }
        }
        if (order == 2 || order == 3)
        {
            if (dimstf[0] != markersdim)
            {
                sprintf(buf, "dimensions of '%s_tf' in prj %d chr %d and alleles array don't match (%lli %lli)\n", DATAARRAY, prj, n, dimstf[0], markersdim);
                printout(buf);
                return -1;
            }
        }
        H5Gclose(chromosome);
    }
    n--;
    if (test == 1)return 1;
    if (nchrinfo == 0)
    {
        chrinfo = (chrinfo_type *) malloc(n * sizeof (chrinfo_type));
        if (chrinfo == NULL)error("malloc error", __LINE__);
    }
    else
    {
        chrinfo = realloc(chrinfo, (nchrinfo + n) * sizeof (chrinfo_type));
        if (chrinfo == NULL)error("realloc error", __LINE__);
    }
    for (i = 0; i < n; i++)
    {
        sprintf(buf, "chr%d", i + 1);
        chromosome = H5Gopen(project, buf, H5P_DEFAULT);
        if (read_str_attr(chromosome, "name", buf) < 0)
        {
            sprintf(chrinfo[nchrinfo + i].name, "chr%d", i + 1);
        }
        else
        {
            sprintf(chrinfo[nchrinfo + i].name, "%s", buf);
        }
        if (order == 1 || order == 3)
        {
            sprintf(buf, "%s_pf", DATAARRAY);
            dataset = H5Dopen(chromosome, buf, H5P_DEFAULT);
            filespace = H5Dget_space(dataset);
            ret = H5Sget_simple_extent_dims(filespace, dimspf, NULL);
            H5Sclose(filespace);
            H5Dclose(dataset);
            dimss[0] = dimspf[1];
            dimss[1] = dimspf[0];
        }
        else
        {
            sprintf(buf, "%s_tf", DATAARRAY);
            dataset = H5Dopen(chromosome, buf, H5P_DEFAULT);
            filespace = H5Dget_space(dataset);
            ret = H5Sget_simple_extent_dims(filespace, dimstf, NULL);
            H5Sclose(filespace);
            H5Dclose(dataset);
            dimss[0] = dimstf[0];
            dimss[1] = dimstf[1];
        }
        chrinfo[nchrinfo + i].dim0 = dimss[0]; //taxa dim
        chrinfo[nchrinfo + i].dim1 = dimss[1]; //positions dim
        chrinfo[nchrinfo + i].mnt = mnt;
        chrinfo[nchrinfo + i].prj = prj;
        chrinfo[nchrinfo + i].num = i + 1;
        H5Gclose(chromosome);
    }
    nchrinfo += n;
    return order;
}

void set_pass()
{
    unsigned long seed[2];
    char salt[] = "$1$........";
    const char *const seedchars = "./0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
    char *password;
    int i;
    FILE *in, *out;

    fprintf(outs, "Password not set - please enter password\n");
    if (!fgets(prtbuf, 100, ins))error("EOF", __LINE__);

    seed[0] = time(NULL);
    srand(seed[0]);
    seed[1] = rand();

    /* Turn it into printable characters from `seedchars'. */
    for (i = 0; i < 8; i++)
        salt[3 + i] = seedchars[(seed[i / 5] >> (i % 5)*6) & 0x3f];

    /* Read in the user's password and encrypt it. */
    password = crypt(prtbuf, salt);
    strcpy(pass, password);
    if ((in = fopen("config.txt", "r")) == NULL)
    {
        fprintf(stderr, "Cannot open config.txt (2)\n");
        exit(1);
    }
    if ((out = fopen("config.txt.new", "w")) == NULL)
    {
        fprintf(stderr, "Cannot open config.txt.new\n");
        exit(1);
    }
    fprintf(out, "#current password\n");
    fprintf(out, "%s\n", pass);
    i = 0;
    while (fgets(prtbuf, 4000, in) != NULL)
    {
        i++;
        if (i > 2)fprintf(out, "%s", prtbuf);
    }
    fclose(in);
    fclose(out);
    if (system("mv -f config.txt.new config.txt") == -1)error("system() call failed", __LINE__);
}

int check_pass(char *testpass)
{
    char *ctestpass;

    ctestpass = crypt(testpass, pass);
    if (strcmp(ctestpass, pass) == 0)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

void set_user_pass(int uid, char *passwd)
{
    unsigned long seed[2];
    char salt[] = "$1$........";
    const char *const seedchars = "./0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
    char *password;
    int i;

    seed[0] = time(NULL);
    srand(seed[0]);
    seed[1] = rand();

    /* Turn it into printable characters from `seedchars'. */
    for (i = 0; i < 8; i++)
        salt[3 + i] = seedchars[(seed[i / 5] >> (i % 5)*6) & 0x3f];

    /* Read in the user's password and encrypt it. */
    password = crypt(passwd, salt);
    strcpy(users[uid].password, password);
}

int check_user_pass(char *name, char *passwd)
{
    int i, n = -1;
    char *passhash;

    if (nusers == 0)return -1;
    for (i = 0; i < nusers; i++)
    {
        //sprintf(prtbuf, "%d %s %s\n", i, users[i].username, users[i].password);
        //printout(prtbuf);
        if (strcmp(users[i].username, name) == 0)n = i;
    }
    if (n == -1)return -1;
    //printout("user OK");
    passhash = crypt(passwd, users[n].password);
    if (strcmp(passhash, users[n].password) == 0)return n;
    return -1;
}

int get_user(char *name)
{
    int i, n = -1;

    if (nusers == 0)return -1;
    for (i = 0; i < nusers; i++)
    {
        //sprintf(prtbuf, "%d %s %s\n", i, users[i].username, users[i].password);
        //printout(prtbuf);
        if (strcmp(users[i].username, name) == 0)n = i;
    }
    return n;
}

hid_t create_udatatype()
{
    hid_t datatype, strtype;

    userdata_type tmppr;
    datatype = H5Tcreate(H5T_COMPOUND, sizeof tmppr);
    if (datatype < 0)error("Cannot create datatype (create_fdatatype)", __LINE__);
    strtype = H5Tcopy(H5T_C_S1);
    H5Tset_size(strtype, 50);
    if (H5Tinsert(datatype, "username", HOFFSET(userdata_type, username), strtype) < 0)error("Cannot insert datatype offset", __LINE__);
    if (H5Tinsert(datatype, "password", HOFFSET(userdata_type, password), strtype) < 0)error("Cannot insert datatype offset", __LINE__);

    return datatype;
}

hid_t create_fdatatype()
{
    hid_t datatype, strtype;

    files_type tmppr;
    datatype = H5Tcreate(H5T_COMPOUND, sizeof tmppr);
    if (datatype < 0)error("Cannot create datatype (create_fdatatype)", __LINE__);
    strtype = H5Tcopy(H5T_C_S1);
    H5Tset_size(strtype, 256);
    if (H5Tinsert(datatype, "name", HOFFSET(files_type, name), strtype) < 0)error("Cannot insert datatype offset", __LINE__);
    if (H5Tinsert(datatype, "mnt", HOFFSET(files_type, mnt), H5T_NATIVE_INT) < 0)error("Cannot insert datatype offset", __LINE__);

    return datatype;
}

hid_t create_pdatatype()
{
    hid_t datatype, strtype;

    projects_type tmppr;
    datatype = H5Tcreate(H5T_COMPOUND, sizeof tmppr);
    if (datatype < 0)error("Cannot create datatype (create_pdatatype)", __LINE__);
    strtype = H5Tcopy(H5T_C_S1);
    H5Tset_size(strtype, 256);
    if (H5Tinsert(datatype, "name", HOFFSET(projects_type, name), strtype) < 0)error("Cannot insert datatype offset", __LINE__);
    if (H5Tinsert(datatype, "num", HOFFSET(projects_type, num), H5T_NATIVE_INT) < 0)error("Cannot insert datatype offset", __LINE__);
    if (H5Tinsert(datatype, "mnt", HOFFSET(projects_type, mnt), H5T_NATIVE_INT) < 0)error("Cannot insert datatype offset", __LINE__);
    if (H5Tinsert(datatype, "order", HOFFSET(projects_type, order), H5T_NATIVE_INT) < 0)error("Cannot insert datatype offset", __LINE__);
    if (H5Tinsert(datatype, "enc", HOFFSET(projects_type, enc), H5T_NATIVE_INT) < 0)error("Cannot insert datatype offset", __LINE__);
    H5Tset_size(strtype, 50);
    if (H5Tinsert(datatype, "encoding", HOFFSET(projects_type, encoding), strtype) < 0)error("Cannot insert datatype offset", __LINE__);

    return datatype;
}

void create_projects()
{
    hid_t dataset, cparms;
    hid_t datatype, dataspace;
    hsize_t dimsf[1], maxdimsf[1], chunk_dims[1];
    herr_t status;
    projects_type tmpprar[1];

    fprintf(flog, "Creating '/projects'\n");
    fflush(flog);
    dimsf[0] = 1;
    maxdimsf[0] = H5S_UNLIMITED;
    chunk_dims[0] = 1;
    dataspace = H5Screate_simple(1, dimsf, maxdimsf);
    if (dataspace < 0)error("Cannot create dataspace (create_projects)", __LINE__);

    datatype = create_pdatatype();

    cparms = H5Pcreate(H5P_DATASET_CREATE);
    if (cparms < 0)error("Cannot create dataset (00)", __LINE__);
    status = H5Pset_chunk(cparms, 1, chunk_dims);
    if (status < 0)error("Cannot create dataset (0)", __LINE__);

    dataset = H5Dcreate(file, "projects", datatype, dataspace, H5P_DEFAULT, cparms, H5P_DEFAULT);
    if (dataset < 0)error("Cannot create dataset", __LINE__);

    status = H5Dextend(dataset, dimsf);
    if (status < 0)error("Cannot extend dataset", __LINE__);
    projects = (projects_type *) malloc(sizeof (projects_type));
    projects[0].mnt = -1;
    projects[0].num = -1;
    projects[0].order = -1;
    projects[0].enc = -1;
    memset(projects[0].name, '\0', 256);
    tmpprar[0] = projects[0];
    status = H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmpprar);
    if (status < 0)error("Cannot write dataset", __LINE__);
    H5Dclose(dataset);
    nprojects = 0;
    fprintf(flog, "'/projects' done\n");
    fflush(flog);
}

void recreate_projects()
{
    hid_t dataset, cparms;
    hid_t datatype, dataspace;
    hsize_t dimsf[1], maxdimsf[1], chunk_dims[1];
    herr_t status;

    fprintf(flog, "Re-Creating '/projects' %d\n", nprojects);
    fflush(flog);
    dimsf[0] = nprojects;
    maxdimsf[0] = H5S_UNLIMITED;
    chunk_dims[0] = 1;
    dataspace = H5Screate_simple(1, dimsf, maxdimsf);
    if (dataspace < 0)error("Cannot create dataspace (create_projects)", __LINE__);

    datatype = create_pdatatype();

    cparms = H5Pcreate(H5P_DATASET_CREATE);
    if (cparms < 0)error("Cannot create dataset (00)", __LINE__);
    status = H5Pset_chunk(cparms, 1, chunk_dims);
    if (status < 0)error("Cannot create dataset (0)", __LINE__);

    dataset = H5Dcreate(file, "projects", datatype, dataspace, H5P_DEFAULT, cparms, H5P_DEFAULT);
    if (dataset < 0)error("Cannot create dataset", __LINE__);

    status = H5Dextend(dataset, dimsf);
    if (status < 0)error("Cannot extend dataset", __LINE__);
    status = H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, projects);
    if (status < 0)error("Cannot write dataset", __LINE__);
    H5Dclose(dataset);
    fprintf(flog, "'/projects' done\n");
    fflush(flog);
}

void create_files()
{
    hid_t dataset, cparms;
    hid_t datatype, dataspace;
    hsize_t dimsf[1], maxdimsf[1], chunk_dims[1];
    herr_t status;

    fprintf(flog, "Creating '/files'\n");
    fflush(flog);
    dimsf[0] = 1;
    maxdimsf[0] = H5S_UNLIMITED;
    chunk_dims[0] = 1;
    dataspace = H5Screate_simple(1, dimsf, maxdimsf);
    if (dataspace < 0)error("Cannot create dataspace (create_files)", __LINE__);
    files_type tmpprar[1];

    datatype = create_fdatatype();

    cparms = H5Pcreate(H5P_DATASET_CREATE);
    status = H5Pset_chunk(cparms, 1, chunk_dims);

    dataset = H5Dcreate(file, "/files", datatype, dataspace, H5P_DEFAULT, cparms, H5P_DEFAULT);
    if (dataset < 0)error("Cannot create dataset", __LINE__);
    status = H5Dextend(dataset, dimsf);
    if (status < 0)error("Cannot extend dataset", __LINE__);
    files = (files_type *) malloc(sizeof (files_type));
    filesh = (hid_t *) malloc(sizeof (hid_t));
    files[0].mnt = -1;
    memset(files[0].name, '\0', 256);
    tmpprar[0] = files[0];
    status = H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmpprar);
    if (status < 0)error("Cannot write dataset", __LINE__);
    H5Dclose(dataset);
    nfiles = 0;
    fprintf(flog, "'/files'  done\n");
    fflush(flog);
}

void recreate_files()
{
    hid_t dataset, cparms;
    hid_t datatype, dataspace;
    hsize_t dimsf[1], maxdimsf[1], chunk_dims[1];
    herr_t status;

    fprintf(flog, "Re-creating '/files'  %d\n", nfiles);
    fflush(flog);
    dimsf[0] = nfiles;
    maxdimsf[0] = H5S_UNLIMITED;
    chunk_dims[0] = 1;
    dataspace = H5Screate_simple(1, dimsf, maxdimsf);
    if (dataspace < 0)error("Cannot create dataspace (create_files)", __LINE__);

    datatype = create_fdatatype();

    cparms = H5Pcreate(H5P_DATASET_CREATE);
    status = H5Pset_chunk(cparms, 1, chunk_dims);

    dataset = H5Dcreate(file, "/files", datatype, dataspace, H5P_DEFAULT, cparms, H5P_DEFAULT);
    if (dataset < 0)error("Cannot create dataset", __LINE__);
    status = H5Dextend(dataset, dimsf);
    if (status < 0)error("Cannot extend dataset", __LINE__);
    status = H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, files);
    if (status < 0)error("Cannot write dataset", __LINE__);
    H5Dclose(dataset);
    fprintf(flog, "'/files'  done\n");
    fflush(flog);
}

void create_users()
{
    hid_t dataset, cparms;
    hid_t datatype, dataspace;
    hsize_t dimsf[1], maxdimsf[1], chunk_dims[1];
    herr_t status;

    fprintf(flog, "Creating '/users'\n");
    fflush(flog);
    dimsf[0] = 1;
    maxdimsf[0] = H5S_UNLIMITED;
    chunk_dims[0] = 1;
    dataspace = H5Screate_simple(1, dimsf, maxdimsf);
    if (dataspace < 0)error("Cannot create dataspace (create_users)", __LINE__);
    userdata_type tmpprar[1];

    datatype = create_udatatype();

    cparms = H5Pcreate(H5P_DATASET_CREATE);
    status = H5Pset_chunk(cparms, 1, chunk_dims);

    dataset = H5Dcreate(file, "/users", datatype, dataspace, H5P_DEFAULT, cparms, H5P_DEFAULT);
    if (dataset < 0)error("Cannot create dataset", __LINE__);
    status = H5Dextend(dataset, dimsf);
    if (status < 0)error("Cannot extend dataset", __LINE__);
    users = (userdata_type *) malloc(sizeof (userdata_type));
    accesstable = (char **) malloc(sizeof (char *));
    memset(users[0].username, '\0', 50);
    accesstable[0] = NULL;
    tmpprar[0] = users[0];
    status = H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmpprar);
    if (status < 0)error("Cannot write dataset", __LINE__);
    nusers = 1;
    if (write_access_table() == -1)
    {
        error("Cannot write accesstable dataset", __LINE__);
    }
    H5Dclose(dataset);
    nusers = 0;
    fprintf(flog, "'/users'  done\n");
    fflush(flog);
}

void recreate_users()
{
    hid_t dataset, cparms;
    hid_t datatype, dataspace;
    hsize_t dimsf[1], maxdimsf[1], chunk_dims[1];
    herr_t status;

    fprintf(flog, "Re-Creating '/users' %d\n", nusers);
    fflush(flog);
    dimsf[0] = nusers;
    maxdimsf[0] = H5S_UNLIMITED;
    chunk_dims[0] = 1;
    dataspace = H5Screate_simple(1, dimsf, maxdimsf);
    if (dataspace < 0)error("Cannot create dataspace (create_users)", __LINE__);

    datatype = create_udatatype();

    cparms = H5Pcreate(H5P_DATASET_CREATE);
    status = H5Pset_chunk(cparms, 1, chunk_dims);

    dataset = H5Dcreate(file, "/users", datatype, dataspace, H5P_DEFAULT, cparms, H5P_DEFAULT);
    if (dataset < 0)error("Cannot create dataset", __LINE__);
    status = H5Dextend(dataset, dimsf);
    if (status < 0)error("Cannot extend dataset", __LINE__);
    status = H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, users);
    if (status < 0)error("Cannot write dataset", __LINE__);
    nusers = 1;
    if (write_access_table() == -1)
    {
        error("Cannot write accesstable dataset", __LINE__);
    }
    H5Dclose(dataset);
    fprintf(flog, "'/users'  done\n");
    fflush(flog);
}

void read_root_props(char c)
{
    hid_t dataset, cparms, datatype;
    hid_t filespace;
    hid_t memspace;
    hsize_t dims[1];
    hsize_t chunk_dims[1];
    herr_t status, status_n;

    fprintf(flog, "Reading root file properties\n");
    fflush(flog);
    int rank, rank_chunk;
    if (c == 'F')
    {
        dataset = H5Dopen(file, "/files", H5P_DEFAULT);
    }
    else if (c == 'P')
    {
        dataset = H5Dopen(file, "/projects", H5P_DEFAULT);
    }
    else // 'U'
    {
        dataset = H5Dopen(file, "/users", H5P_DEFAULT);
    }
    if (dataset < 0)error("Cannot open dataset", __LINE__);

    filespace = H5Dget_space(dataset); /* Get filespace handle first. */
    rank = H5Sget_simple_extent_ndims(filespace);
    status_n = H5Sget_simple_extent_dims(filespace, dims, NULL);
    cparms = H5Dget_create_plist(dataset); /* Get properties handle first. */
    rank_chunk = H5Pget_chunk(cparms, 1, chunk_dims);
    memspace = H5Screate_simple(1, dims, NULL);
    fprintf(flog, " = %c %lli\n", c, dims[0]);
    if (c == 'F')
    {
        if (nfiles > 0)
        {
            free(files);
            free(filesh);
        }
        files = (files_type *) malloc(dims[0] * sizeof (files_type));
        filesh = (hid_t *) malloc(dims[0] * sizeof (hid_t));
        datatype = create_fdatatype();
        status = H5Dread(dataset, datatype, memspace, filespace, H5P_DEFAULT, files);
        if (status < 0)error("Cannot read files dataset", __LINE__);
        nfiles = dims[0];
        if (files[0].mnt == -1)nfiles = 0;
    }
    else if (c == 'P')
    {
        if (nprojects > 0)
        {
            free(projects);
        }
        projects = (projects_type *) malloc(dims[0] * sizeof (projects_type));
        datatype = create_pdatatype();
        status = H5Dread(dataset, datatype, memspace, filespace, H5P_DEFAULT, projects);
        if (status < 0)error("Cannot read projects dataset", __LINE__);
        nprojects = dims[0];
        if (projects[0].mnt == -1)nprojects = 0;
    }
    else // 'U'
    {
        if (nusers > 0)
        {
            free(users);
            freestrtable(accesstable, nusers);
        }
        users = (userdata_type *) malloc(dims[0] * sizeof (userdata_type));
        datatype = create_udatatype();
        status = H5Dread(dataset, datatype, memspace, filespace, H5P_DEFAULT, users);
        if (status < 0)error("Cannot read users dataset", __LINE__);
        nusers = dims[0];
        read_access_table();
        if (users[0].username[0] == '\0')nusers = 0;
    }

    H5Dclose(dataset);
}

void write_root_props(char c, int shrink)
{
    hid_t dataset, cparms, datatype;
    hid_t filespace;
    hid_t memspace;
    hsize_t dims[1], dimsf[1];
    hsize_t chunk_dims[1];
    herr_t status, status_n;

    fprintf(flog, "Writing root file properties\n");
    fflush(flog);
    int rank, rank_chunk;
    if (c == 'F')
    {
        if (shrink == 1)
        {
            status = H5Ldelete(file, "/files", H5P_DEFAULT);
            recreate_files();
            return;
        }
        dataset = H5Dopen(file, "/files", H5P_DEFAULT);
    }
    else if (c == 'P')
    {
        if (shrink == 1)
        {
            status = H5Ldelete(file, "/projects", H5P_DEFAULT);
            recreate_projects();
            return;
        }
        dataset = H5Dopen(file, "/projects", H5P_DEFAULT);
    }
    else if (c == 'U')
    {
        if (shrink == 1)
        {
            status = H5Ldelete(file, "/users", H5P_DEFAULT);
            status = H5Ldelete(file, "/access", H5P_DEFAULT);
            recreate_users();
            return;
        }
        dataset = H5Dopen(file, "/users", H5P_DEFAULT);
    }
    else // 'A'
    {
        status = H5Ldelete(file, "/access", H5P_DEFAULT);
        if (write_access_table() == -1)error("Cannot write users access dataset", __LINE__);
        return;
    }
    if (dataset < 0)error("Cannot open dataset", __LINE__);

    filespace = H5Dget_space(dataset); /* Get filespace handle first. */
    rank = H5Sget_simple_extent_ndims(filespace);
    status_n = H5Sget_simple_extent_dims(filespace, dims, NULL);
    if (c == 'F')if (dims[0] > nfiles)dims[0] = (hsize_t) nfiles;
    if (c == 'P')if (dims[0] > nprojects)dims[0] = (hsize_t) nprojects;
    if (c == 'U')if (dims[0] > nusers)dims[0] = (hsize_t) nusers;
    cparms = H5Dget_create_plist(dataset); /* Get properties handle first. */
    rank_chunk = H5Pget_chunk(cparms, 1, chunk_dims);
    memspace = H5Screate_simple(1, dims, NULL);
    fprintf(flog, " = %c %lli\n", c, dims[0]);
    fflush(flog);
    if (c == 'F')
    {
        datatype = create_fdatatype();
        if (nfiles > dims[0] && (files[0].mnt != -1))
        {
            //increase the size of dataset in file
            cparms = H5Pcreate(H5P_DATASET_CREATE);
            status = H5Pset_chunk(cparms, 1, chunk_dims);
            dimsf[0] = nfiles;
            status = H5Dextend(dataset, dimsf);
            if (status < 0)error("Cannot extend dataset", __LINE__);
        }
        status = H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, files);
        if (status < 0)error("Cannot write files dataset", __LINE__);
    }
    else if (c == 'P')
    {
        if (nprojects > dims[0] && (projects[0].mnt != -1))
        {
            //increase the size of dataset in file
            cparms = H5Pcreate(H5P_DATASET_CREATE);
            status = H5Pset_chunk(cparms, 1, chunk_dims);
            dimsf[0] = nprojects;
            status = H5Dextend(dataset, dimsf);
            if (status < 0)error("Cannot extend dataset", __LINE__);
        }
        datatype = create_pdatatype();
        status = H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, projects);
        if (status < 0)error("Cannot write projects dataset", __LINE__);
    }
    else //'U'
    {
        if (nusers > dims[0] && (users[0].username[0] != '\0'))
        {
            //increase the size of dataset in file
            cparms = H5Pcreate(H5P_DATASET_CREATE);
            status = H5Pset_chunk(cparms, 1, chunk_dims);
            dimsf[0] = nusers;
            status = H5Dextend(dataset, dimsf);
            if (status < 0)error("Cannot extend dataset", __LINE__);
        }
        datatype = create_udatatype();
        status = H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, users);
        if (status < 0)error("Cannot write users dataset", __LINE__);
    }

    H5Dclose(dataset);
}

int main(int argc, char** argv)
{
    char command[50];
    FILE *in, *fptr;
    time_t ttt;
    int pid, iii;

    nfiles = 0;
    nprojects = 0;
    nchrinfo = 0;
    comm_pid = -1;

    strcpy(hdf5conf, INSTALLDIR);
    if (argc == 2)
    {
        strcpy(hdf5conf, argv[1]);
    }

    //open log file
    sprintf(prtbuf, "%s/cbsuhdf5.log", hdf5conf);
    if ((flog = fopen(prtbuf, "a")) == NULL)
    {
        fprintf(stderr, "ERROR: Cannot open hdf5.log\n");
        exit(1);
    }
    pid = getpid();
    //create pid file
    sprintf(prtbuf, "%s/cbsuhdf5.pid", hdf5conf);
    if ((fptr = fopen(prtbuf, "w")) == NULL)
    {
        fprintf(stderr, "ERROR: Cannot open cbsuhdf5.pid\n");
        exit(1);
    }
    fprintf(fptr, "%d\n", pid);
    fclose(fptr);

    //read config file first
    sprintf(prtbuf, "%s/config.txt", hdf5conf);
    delayed = 1;
    outsbuffer = NULL;
    if ((in = fopen(prtbuf, "r")) != NULL)
    {
        if (!fgets(prtbuf, 4000, in))error("EOF", __LINE__);
        if (!fgets(pass, 100, in))error("EOF", __LINE__);
        chomp(pass);
        if (!fgets(prtbuf, 4000, in))error("EOF", __LINE__);
        if (!fgets(hdf5dir, 500, in))error("EOF", __LINE__);
        chomp(hdf5dir);
        ttt = time(NULL);
        fprintf(flog, "HDF5 CBSU Database Server %s\n", asctime(localtime(&ttt)));
        fflush(flog);
        if (!fgets(prtbuf, 4000, in))error("EOF", __LINE__);
        if (!fgets(hdf5root, 500, in))error("EOF", __LINE__);
        chomp(hdf5root);
        sprintf(prtbuf, "%s%s", hdf5dir, hdf5root);
        strcpy(hdf5root, prtbuf);
        if (!fgets(prtbuf, 4000, in))error("EOF", __LINE__);
        if (!fgets(tmpdir, 500, in))error("EOF", __LINE__);
        chomp(tmpdir);
        if (!fgets(prtbuf, 4000, in))error("EOF", __LINE__);
        if (!fgets(command, 50, in))error("EOF", __LINE__);
        chomp(command);
        fclose(in);
        //set communication mode
        if (strcmp(command, "std") == 0)
        {
            outs = stdout;
            ins = stdin;
        }
        else if (command[0] == '_')
        {
            //open named pipes here 
            sprintf(prtbuf, "rm -f %s%s_out", hdf5dir, &command[1]);
            if (system(prtbuf) == -1)error("system() call failed", __LINE__);
            sprintf(prtbuf, "rm -f %s%s_in", hdf5dir, &command[1]);
            if (system(prtbuf) == -1)error("system() call failed", __LINE__);
            sprintf(prtbuf, "%s%s_out", hdf5dir, &command[1]);
            int ret_val = mkfifo(prtbuf, 0666);
            if ((ret_val == -1) && (errno != EEXIST))
            {
                fprintf(flog, "Error creating the named pipe %s\n", prtbuf);
                fprintf(stderr, "Error creating the named pipe %s\n", prtbuf);
                fclose(flog);
                exit(1);
            }
            sprintf(prtbuf, "%s%s_in", hdf5dir, &command[1]);
            ret_val = mkfifo(prtbuf, 0666);
            if ((ret_val == -1) && (errno != EEXIST))
            {
                fprintf(flog, "Error creating the named pipe %s\n", prtbuf);
                fprintf(stderr, "Error creating the named pipe %s\n", prtbuf);
                fclose(flog);
                exit(1);
            }
            //creating child comm program
            comm_pid = fork();
            if (comm_pid < 0)
            {
                error("fork failed\n", __LINE__);
            }
            if (comm_pid == 0)
            {
                sprintf(prtbuf, "%s%s_out", hdf5dir, &command[1]);
                fprintf(flog, "COMM: Opening %s pipe for reading\n", prtbuf);
                fflush(flog);
                if ((outs = fopen(prtbuf, "r")) == NULL)
                {
                    fprintf(stderr, "COMM: Cannot open output pipe descriptor %s\n", prtbuf);
                    exit(1);
                }
                sprintf(prtbuf, "%s%s_in", hdf5dir, &command[1]);
                fprintf(flog, "COMM Opening %s pipe for writing\n", prtbuf);
                fflush(flog);
                if ((ins = fopen(prtbuf, "w")) == NULL)
                {
                    fprintf(stderr, "COMM: Cannot open input pipe descriptor %s\n", prtbuf);
                    exit(1);
                }
                while (1)
                {

                }
                exit(0);
            }
            delayed = 1;
            //create pid file
            sprintf(prtbuf, "%s/cbsuhdf5.pid1", hdf5conf);
            if ((fptr = fopen(prtbuf, "w")) == NULL)
            {
                fprintf(stderr, "ERROR: Cannot open cbsuhdf5.pid1\n");
                exit(1);
            }
            fprintf(fptr, "%d\n", comm_pid);
            fclose(fptr);
            sprintf(prtbuf, "%s%s_out", hdf5dir, &command[1]);
            fprintf(flog, "Opening %s pipe for writing\n", prtbuf);
            fflush(flog);
            if ((outs = fopen(prtbuf, "w")) == NULL)
            {
                fprintf(stderr, "Cannot open output pipe descriptor %s\n", prtbuf);
                exit(1);
            }
            sprintf(prtbuf, "%s%s_in", hdf5dir, &command[1]);
            fprintf(flog, "Opening %s pipe for reading\n", prtbuf);
            fflush(flog);
            if ((ins = fopen(prtbuf, "r")) == NULL)
            {
                fprintf(stderr, "Cannot open input pipe descriptor %s\n", prtbuf);
                exit(1);
            }
            sprintf(prtbuf, "rm -f %shdf5.lock", hdf5dir);
            iii = system(prtbuf);
        }
        else
        {
            //invalid option
            outs = stdout;
            ins = stdin;
            fprintf(stderr, "INVALID communication option %s, using std instead\n", command);
        }
    }
    else
    {
        fprintf(stderr, "ERROR: Cannot open config.txt\n");
        exit(1);
    }
    //if database file does not exists, initialize it
    struct stat stt;
    int i = stat(hdf5root, &stt);
    int j = 0;
    fflush(flog);
    if (i != 0)
    {
        fprintf(flog, "Database %s does not exist - initializing\n", hdf5root);
        fflush(flog);
        file = H5Fcreate(hdf5root, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        //create storage for files and projects
        create_projects();
        create_files();
        create_users();
        H5Fflush(file, H5F_SCOPE_LOCAL);
    }
    else
    {
        fprintf(flog, "Opening existing database %s\n", hdf5root);
        fflush(flog);
        file = H5Fopen(hdf5root, H5F_ACC_RDWR, H5P_DEFAULT);
        if (file < 0)error("Cannot open root hdf5 file", __LINE__);
        //read in all previously mounted files and remount them
        read_root_props('F');
        read_root_props('P');
        read_root_props('U');
        fprintf(flog, "nfiles=%d nprojects=%d nusers=%d\n", nfiles, nprojects, nusers);
        fflush(flog);
        for (i = 0; i < nfiles; i++)
        {
            //mount file
            if (files[i].mnt != -1)
            {
                fprintf(flog, "Mounting file %s\n", files[i].name);
                fflush(flog);
                sprintf(prtbuf, "%s%s", hdf5dir, files[i].name);
                filesh[i] = H5Fopen(prtbuf, H5P_DEFAULT, H5P_DEFAULT);
                sprintf(command, "/mnt%d", files[i].mnt);
                if (H5Fmount(file, command, filesh[i], H5P_DEFAULT) < 0)
                {
                    sprintf(prtbuf, "Cannot mount file %s\n", files[i].name);
                    //remove file from cache
                    umount1(i);
                    sprintf(prtbuf, "File %s has been removed from configuration\nPlease restart the server\n", files[i].name);
                    error(prtbuf, __LINE__);
                }
                for (j = 0; j < nprojects; j++)
                {
                    if (projects[j].mnt == i)
                    {
                        sprintf(command, "/mnt%d/project%d", i, projects[j].num);
                        hid_t project = H5Gopen(file, command, H5P_DEFAULT);
                        if (project < 0)
                        {
                            sprintf(prtbuf, "Cannot open %s\n", command);
                            delproject(j);
                            sprintf(prtbuf, "Project %d has been removed from configuration\nPlease restart the server\n", j);
                            error(prtbuf, __LINE__);
                        }
                        if (read_chrinfo(project, i, j, 0, 0, projects[j].order) < 0)
                        {
                            sprintf(prtbuf, "Error reading chromosomes from %s\n", command);
                            error(prtbuf, __LINE__);
                        }
                        H5Gclose(project);
                    }
                }
            }
        }
    }

    //check if the password id not empty
    if (strlen(pass) == 0)
    {
        //set password
        fprintf(flog, "Setting admin password\n");
        set_pass();
    }

    //enter command loop, command ends on empty line
    fprintf(flog, "entering command loop\n\n");
    fflush(flog);
    if (delayed == 1)
    {
        if (outsbuffer != NULL)
        {
            free(outsbuffer);
            outsbuffer = NULL;
        }
    }
    while (1)
    {
        errbuf[0] = '\0';
        if (!fgets(command, 50, ins))error("EOF", __LINE__);
        chomp(command);
        fprintf(flog, "Command: %s\n", command);
        ttt = time(NULL);
        fprintf(flog, "%s", asctime(localtime(&ttt)));
        fflush(flog);
        gettimeofday(&c_start0, NULL);
        gettimeofday(&c_start, NULL);
        int iret = -99;
        if (strcmp(command, "QUERY") == 0)
        {
            iret = query();
        }
        else if (strcmp(command, "TABLE") == 0)
        {
            iret = table();
        }
        else if (strcmp(command, "MOUNT") == 0)
        {
            iret = mount();
        }
        else if (strcmp(command, "UMOUNT") == 0)
        {
            iret = umount();
        }
        else if (strcmp(command, "INDEX") == 0)
        {
            iret = sort();
        }
        else if (strcmp(command, "PINFO") == 0)
        {
            iret = pinfo();
        }
        else if (strcmp(command, "PLIST") == 0)
        {
            iret = plist();
        }
        else if (strcmp(command, "FLIST") == 0)
        {
            iret = flist();
        }
        else if (strcmp(command, "AFLIST") == 0)
        {
            iret = aflist();
        }
        else if (strcmp(command, "FINFO") == 0)
        {
            iret = finfo();
        }
        else if (strcmp(command, "LISTALL") == 0)
        {
            iret = cache();
        }
        else if (strcmp(command, "LOGIN") == 0)
        {
            iret = login();
        }
        else if (strcmp(command, "USERADD") == 0)
        {
            iret = useradd();
        }
        else if (strcmp(command, "USERDEL") == 0)
        {
            iret = userdel();
        }
        else if (strcmp(command, "USERPASS") == 0)
        {
            iret = userpass();
        }
        else if (strcmp(command, "USERACC") == 0)
        {
            iret = useracc();
        }
        else if (strcmp(command, "ULIST") == 0)
        {
            iret = ulist();
        }
        else if (strcmp(command, "QUIT") == 0)
        {
            iret = quit();
            if (iret == 1)
            {
                printout("Command COMPLETED\n");
                break;
            }
        }
        else
        {
            printout("INVALID COMMAND\n");
            iret = 0;
        }
        if (errbuf[0] != '\0')printout(errbuf);
        gettimeofday(&c_end, NULL);
        if (iret == 0 && errbuf[0] != '\0')printout(errbuf);
        if (delayed == 1)printoutbuf();
        int delayed1 = delayed;
        delayed = 0;
        ttt = time(NULL);
        sprintf(prtbuf, "END %s", asctime(localtime(&ttt)));
        printout(prtbuf);
        if (iret == -1)
        {
            printout("Incomplete command sequence\n");
        }
        else if (iret == 0)
        {
            printout("Command execution error\n");
        }
        else
        {
            printout("Command execution successful\n");
        }
        double etime = ((double) c_end.tv_sec + ((double) c_end.tv_usec) / 1000000)-((double) c_start.tv_sec + ((double) c_start.tv_usec) / 1000000);
        sprintf(prtbuf, "execution %10.3lf seconds\n", etime);
        printout(prtbuf);
        etime = ((double) c_end.tv_sec + ((double) c_end.tv_usec) / 1000000)-((double) c_start0.tv_sec + ((double) c_start0.tv_usec) / 1000000);
        sprintf(prtbuf, "total     %10.3lf seconds\n", etime);
        printout(prtbuf);
        printout("Command COMPLETED\n");
        delayed = delayed1;
    }
    ttt = time(NULL);
    fprintf(flog, "Exiting %s\n", asctime(localtime(&ttt)));
    fclose(flog);
    return (EXIT_SUCCESS);
}

void create_attrib_str(hid_t handle, char *name, char* value)
{
    hid_t aid3, atype, attr3;
    herr_t ret;

    aid3 = H5Screate(H5S_SCALAR);
    atype = H5Tcopy(H5T_C_S1);
    H5Tset_size(atype, strlen(value));
    attr3 = H5Acreate(handle, name, atype, aid3, H5P_DEFAULT, H5P_DEFAULT);
    ret = H5Awrite(attr3, atype, value);
    ret = H5Sclose(aid3);
    ret = H5Aclose(attr3);

}

int read_int(hsize_t ioff, hsize_t icount, hid_t dataset, hid_t datatype, hid_t dataspace, int **valout)
{
    hsize_t memdim[2], count[2], offset[2], memoffset[2];
    herr_t ret;
    hid_t memspace;
    int *val;

    offset[0] = ioff;
    count[0] = icount;
    memoffset[0] = 0;
    memdim[0] = icount;

    ret = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);
    if (ret < 0)
    {
        sprintf(prtbuf, "Cannot select hyperslab in database\n");
        printout(prtbuf);
        return -1;
    }
    memspace = H5Screate_simple(1, memdim, NULL);
    ret = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, memoffset, NULL, count, NULL);
    if (ret < 0)
    {
        sprintf(prtbuf, "Cannot select hyperslab in memory\n");
        printout(prtbuf);
        return -1;
    }

    val = malloc(icount * sizeof (int));
    if (val == NULL)
    {
        sprintf(prtbuf, "Cannot allocate memory (%llu)\n", sizeof (int) *icount);
        printout(prtbuf);
        return -1;
    }

    ret = H5Dread(dataset, datatype, memspace, dataspace, H5P_DEFAULT, val);
    if (ret < 0)
    {
        sprintf(prtbuf, "Cannot read data\n");
        printout(prtbuf);
        H5Sclose(memspace);
        free(val);
        return -1;
    }

    *valout = val;

    H5Sclose(memspace);

    return 1;

}

int read_str(hsize_t ioff, hsize_t icount, hid_t dataset, hid_t datatype, hid_t dataspace, char ***valout)
{
    hsize_t memdim[2], count[2], offset[2], memoffset[2];
    herr_t ret;
    hid_t memspace, xfer_pid;
    char **val;

    offset[0] = ioff;
    count[0] = icount;
    memoffset[0] = 0;
    memdim[0] = icount;

    ret = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);
    if (ret < 0)
    {
        sprintf(prtbuf, "Cannot select hyperslab in database\n");
        printout(prtbuf);
        return -1;
    }
    memspace = H5Screate_simple(1, memdim, NULL);
    ret = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, memoffset, NULL, count, NULL);
    if (ret < 0)
    {
        sprintf(prtbuf, "Cannot select hyperslab in memory\n");
        printout(prtbuf);
        return -1;
    }
    xfer_pid = H5Pcreate(H5P_DATASET_XFER);
    if (xfer_pid < 0)
    {
        sprintf(prtbuf, "Cannot get xfrer_pid\n");
        printout(prtbuf);
        return -1;
    }

    val = malloc(icount * sizeof (char *));
    if (val == NULL)
    {
        sprintf(prtbuf, "Cannot allocate memory (%llu)\n", sizeof (char *) *icount);
        printout(prtbuf);
        return -1;
    }

    ret = H5Dread(dataset, datatype, memspace, dataspace, xfer_pid, val);
    if (ret < 0)
    {
        sprintf(prtbuf, "Cannot read data\n");
        printout(prtbuf);
        H5Sclose(memspace);
        free(val);
        return -1;
    }

    H5Sclose(memspace);
    H5Pclose(xfer_pid);

    *valout = val;

    return 1;
}

int readone(hsize_t ioff, hid_t dataset, hid_t datatype, hid_t dataspace, char *val)
{
    hsize_t memdim[2], count[2], offset[2], memoffset[2];
    herr_t ret;
    hid_t memspace;

    offset[0] = ioff;
    count[0] = 1;
    memoffset[0] = 0;
    memdim[0] = 1;

    ret = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);
    if (ret < 0)
    {
        sprintf(prtbuf, "Cannot select hyperslab in database\n");
        printout(prtbuf);
        return -1;
    }
    memspace = H5Screate_simple(1, memdim, NULL);
    ret = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, memoffset, NULL, count, NULL);
    if (ret < 0)
    {
        sprintf(prtbuf, "Cannot select hyperslab in memory\n");
        printout(prtbuf);
        return -1;
    }

    ret = H5Dread(dataset, datatype, memspace, dataspace, H5P_DEFAULT, val);
    if (ret < 0)
    {
        sprintf(prtbuf, "Cannot read data\n");
        printout(prtbuf);
        H5Sclose(memspace);
        return -1;
    }

    H5Sclose(memspace);

    return 1;

}

int readone_str(hsize_t ioff, hid_t dataset, hid_t datatype, hid_t dataspace, char *val)
{
    hsize_t memdim[2], count[2], offset[2], memoffset[2];
    herr_t ret;
    hid_t memspace, xfer_pid;
    char *valptr[2];

    offset[0] = ioff;
    count[0] = 1;
    memoffset[0] = 0;
    memdim[0] = 1;

    ret = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);
    if (ret < 0)
    {
        sprintf(prtbuf, "Cannot select hyperslab in database\n");
        printout(prtbuf);
        return -1;
    }
    memspace = H5Screate_simple(1, memdim, NULL);
    ret = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, memoffset, NULL, count, NULL);
    if (ret < 0)
    {
        sprintf(prtbuf, "Cannot select hyperslab in memory\n");
        printout(prtbuf);
        return -1;
    }
    xfer_pid = H5Pcreate(H5P_DATASET_XFER);
    if (xfer_pid < 0)
    {
        sprintf(prtbuf, "Cannot get xfrer_pid\n");
        printout(prtbuf);
        return -1;
    }

    ret = H5Dread(dataset, datatype, memspace, dataspace, xfer_pid, valptr);
    if (ret < 0)
    {
        sprintf(prtbuf, "Cannot read data\n");
        printout(prtbuf);
        H5Sclose(memspace);
        return -1;
    }

    H5Sclose(memspace);
    H5Pclose(xfer_pid);

    strcpy(val, valptr[0]);

    free(valptr[0]);

    return 1;
}

int binsearch_str(hid_t obj, char *tblname, char *value, long unsigned int *pp)
{
    hid_t dataspace, dataset, datatype;
    long unsigned int i1, i2, ii, kk;
    int status_n, nn;
    char v1[MAXSTRLEN], v2[MAXSTRLEN], vv[MAXSTRLEN];
    hsize_t dims_out[2];
    char indx_i[256], indx_v[256];
    unsigned int jj;

    sprintf(indx_i, "%s_si", tblname);
    if (!H5Lexists(obj, indx_i, H5P_DEFAULT))
    {
        sprintf(prtbuf, "Table '%s' is not indexed, please run INDEX command on this file\n", tblname);
        printout(prtbuf);
        return -1;
    }
    sprintf(indx_v, "%s_sv", tblname);
    if (!H5Lexists(obj, indx_v, H5P_DEFAULT))
    {
        sprintf(prtbuf, "Table '%s' is not indexed, please run INDEX command on this file\n", tblname);
        printout(prtbuf);
        return -1;
    }
    //open dataset
    dataset = H5Dopen(obj, indx_v, H5P_DEFAULT);
    if (dataset < 0)
    {
        sprintf(prtbuf, "Cannot open dataset '%s'\n", indx_v);
        printout(prtbuf);
        return -1;
    }
    dataspace = H5Dget_space(dataset);
    if (dataspace < 0)
    {
        sprintf(prtbuf, "Cannot open dataspace in dataset '%s'\n", indx_v);
        printout(prtbuf);
        H5Dclose(dataset);
        return -1;
    }

    datatype = H5Dget_type(dataset);

    status_n = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);

    i1 = 0;
    if (readone_str(i1, dataset, datatype, dataspace, v1) == -1)
    {
        sprintf(prtbuf, "Cannot read element %lu in dataset '%s'\n", i1, tblname);
        printout(prtbuf);
        H5Dclose(dataset);
        H5Sclose(dataspace);
        return -1;
    }
    i2 = dims_out[0] - 1;
    if (readone_str(i2, dataset, datatype, dataspace, v2) == -1)
    {
        sprintf(prtbuf, "Cannot read element %lu in dataset '%s'\n", i2, tblname);
        printout(prtbuf);
        H5Dclose(dataset);
        H5Sclose(dataspace);
        return -1;
    }
    int search = 1;
    if (strcmp(v1, value) == 0)
    {
        kk = i1;
        search = 0;
    }
    if (strcmp(v2, value) == 0)
    {
        kk = i2;
        search = 0;
    }
    nn = 0;
    while (search)
    {
        nn++;
        //fprintf(flog, "'%s' == %li '%s' %li '%s' iter %d\n", value, i1, v1, i2, v2, nn);
        //fflush(flog);
        if (i1 + 1 == i2 || nn > ITERLIMIT)
        {
            sprintf(prtbuf, "Cannot find value %s in dataset '%s' in %d iterations\n", value, tblname, nn - 1);
            printout(prtbuf);
            H5Dclose(dataset);
            H5Sclose(dataspace);
            return -1;
        }
        ii = (i1 + i2) / 2;
        if (readone_str(ii, dataset, datatype, dataspace, vv) == -1)
        {
            sprintf(prtbuf, "Cannot read element %lu in dataset '%s'\n", ii, indx_v);
            printout(prtbuf);
            H5Dclose(dataset);
            H5Sclose(dataspace);
            return -1;
        }
        //fprintf(flog, "'%s' ==> %li '%s' ", value, ii, vv);
        if (strcmp(vv, value) == 0)
        {
            kk = ii;
            break;
        }
        if (strcmp(value, vv) > 0)
        {
            //fprintf(flog, " RIGHT\n");
            //fflush(flog);
            i1 = ii;
            strcpy(v1, vv);
        }
        else
        {
            //fprintf(flog, " LEFT\n");
            //fflush(flog);
            i2 = ii;
            strcpy(v2, vv);
        }
    }

    H5Dclose(dataset);
    H5Sclose(dataspace);

    dataset = H5Dopen(obj, indx_i, H5P_DEFAULT);
    if (dataset < 0)
    {
        sprintf(prtbuf, "Cannot open dataset '%s'\n", indx_i);
        printout(prtbuf);
        return -1;
    }
    dataspace = H5Dget_space(dataset);
    if (dataspace < 0)
    {
        sprintf(prtbuf, "Cannot open dataspace in dataset '%s'\n", indx_i);
        printout(prtbuf);
        H5Dclose(dataset);
        return -1;
    }

    datatype = H5Dget_type(dataset);

    if (readone(kk, dataset, datatype, dataspace, (char *) &jj) == -1)
    {
        sprintf(prtbuf, "Cannot read element %lu in dataset '%s'\n", kk, indx_i);
        printout(prtbuf);
        H5Dclose(dataset);
        H5Sclose(dataspace);
        return -1;
    }

    *pp = (long unsigned int) jj;
    //fprintf(flog, "== binsearch_str INDEX %u iter %d\n", jj, nn);
    //fflush(flog);
    H5Dclose(dataset);
    H5Sclose(dataspace);
    return 1;
}

int binsearch_int(hid_t obj, char *tblname, char *value, long unsigned int *pp, int round)
{
    hid_t dataspace, dataset, datatype;
    long unsigned int i1, i2, ii;
    int v1, v2, vv, status_n, val = 0, nn;
    hsize_t dims_out[2];

    sscanf(value, "%d", &val);
    //open dataset
    dataset = H5Dopen(obj, tblname, H5P_DEFAULT);
    if (dataset < 0)
    {
        sprintf(prtbuf, "Cannot open dataset '%s'\n", tblname);
        printout(prtbuf);
        return -1;
    }
    dataspace = H5Dget_space(dataset);
    if (dataspace < 0)
    {
        sprintf(prtbuf, "Cannot open dataspace in dataset '%s'\n", tblname);
        printout(prtbuf);
        H5Dclose(dataset);
        return -1;
    }
    datatype = H5Tcopy(H5T_NATIVE_INT);

    status_n = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);

    i1 = 0;
    if (readone(i1, dataset, datatype, dataspace, (char *) &v1) == -1)
    {
        sprintf(prtbuf, "Cannot read element %lu in dataset '%s'\n", i1, tblname);
        printout(prtbuf);
        H5Dclose(dataset);
        H5Sclose(dataspace);
        return -1;
    }
    i2 = dims_out[0] - 1;
    if (readone(i2, dataset, datatype, dataspace, (char *) &v2) == -1)
    {
        sprintf(prtbuf, "Cannot read element %lu in dataset '%s'\n", i2, tblname);
        printout(prtbuf);
        H5Dclose(dataset);
        H5Sclose(dataspace);
        return -1;
    }
    if (v1 == val)
    {
        *pp = i1;
        return 1;
    }
    if (v2 == val)
    {
        *pp = i2;
        return 1;
    }
    nn = 0;
    while (1)
    {
        nn++;
        //fprintf(flog, "%d == %li %d %li %d iter %d\n", val, i1, v1, i2, v2, nn);
        //fflush(flog);
        if (i1 + 1 == i2)
        {
            if (round == -1)
            {
                sprintf(prtbuf, "Cannot find value %s in dataset '%s' in %d iterations\n", value, tblname, nn - 1);
                printout(prtbuf);
                H5Dclose(dataset);
                H5Sclose(dataspace);
                return -1;
            }
            else
            {
                *pp = i2;
                if (round == 0)*pp = i1;
                //fprintf(flog, "== INDEX %lu\n", ii);
                //fflush(flog);
                H5Dclose(dataset);
                H5Sclose(dataspace);
                //fprintf(flog, "== binsearch_int INDEX %lu iter %d ROUNDED %d\n", *pp, nn, round);
                //fflush(flog);
                return 1;
            }
        }
        if (nn > ITERLIMIT)
        {
            sprintf(prtbuf, "Cannot find value %s in dataset '%s' in %d iterations\n", value, tblname, nn - 1);
            printout(prtbuf);
            H5Dclose(dataset);
            H5Sclose(dataspace);
            return -1;
        }
        ii = (i1 + i2) / 2;
        if (readone(ii, dataset, datatype, dataspace, (char *) &vv) == -1)
        {
            sprintf(prtbuf, "Cannot read element %lu in dataset '%s'\n", ii, tblname);
            printout(prtbuf);
            H5Dclose(dataset);
            H5Sclose(dataspace);
            return -1;
        }
        //fprintf(flog, "%d ==> %li %d\n", val, ii, vv);
        //fflush(flog);
        if (vv == val)
        {
            *pp = ii;
            //fprintf(flog, "== INDEX %lu\n", ii);
            //fflush(flog);
            H5Dclose(dataset);
            H5Sclose(dataspace);
            //fprintf(flog, "== binsearch_int INDEX %lu iter %d\n", ii, nn);
            //fflush(flog);
            return 1;
        }
        if (vv < val)
        {
            i1 = ii;
            v1 = vv;
        }
        else
        {
            i2 = ii;
            v2 = vv;
        }
    }


    return -1;
}

int convert_position(char *buf, hid_t chrm, long unsigned int *pp, int i, char *valuetype, char *rangetype)
{
    //find indexes if necessary
    if (strcmp(valuetype, "positions") == 0)
    {
        int round = -1;
        if (strcmp(rangetype, "range") == 0)round = 1 - i;
        return binsearch_int(chrm, "positions", buf, &pp[i], round);
    }
    else if (strcmp(valuetype, "markers") == 0)
    {
        return binsearch_str(chrm, "markers", buf, &pp[i]);
    }
    else if (strcmp(valuetype, "indexes") == 0)
    {
        sscanf(buf, "%li", &pp[i]);
    }
    else
    {
        //unsupported type
        printout("Invalid positions range value type\n");
        return -1;
    }
    return 1;
}

int convert_taxa(char *buf, hid_t prj, long unsigned int *pp, int i, char *valuetype)
{
    //find indexes if necessary
    if (strcmp(valuetype, "taxa") == 0)
    {
        return binsearch_str(prj, "taxa", buf, &pp[i]);
    }
    else if (strcmp(valuetype, "indexes") == 0)
    {
        sscanf(buf, "%li", &pp[i]);
    }
    else
    {
        //unsupported type
        printout("Invalid taxa range value type\n");
        return -1;
    }
    return 1;
}

void printdata(int prj, char *data, int iformat, int ntt, int npp, FILE *outf, int pstride, int tstride, int order)
{
    int ii, jj;

    int enc = projects[prj].enc;
    int npps = npp / pstride;
    int ntts = ntt / tstride;
    int loop1 = npps;
    int loop2 = ntts;
    if (order == 2)
    {
        loop1 = ntts;
        loop2 = npps;
    }

    for (ii = 0; ii < loop2; ii++)
    {
        for (jj = 0; jj < loop1; jj++)
        {
            int offs;
            //if (projects[prj].order == 2)
            //{
            //    offs = ii * ntts + jj;
            //}
            //else
            //{
            //    offs = jj * npps + ii;
            //}
            offs = jj + ii * loop1;
            if (enc == 1)
            {
                char c = data[offs];
                if (c == '\0')c = '.';
                if (iformat == 0)
                {
                    fprintf(outf, "%c", c);
                }
                else
                {
                    if (jj > 0)fprintf(outf, " ");
                    fprintf(outf, "%d", c);
                }
            }
            else
            {
                char c0 = data[offs * enc];
                char c1 = data[offs * enc + 1];
                char c2 = data[offs * enc + 2];
                if (c0 == '\0')c0 = '.';
                if (jj > 0)fprintf(outf, " ");
                if (iformat == 0)
                {
                    fprintf(outf, "%c %d %d", c0, c1, c2);
                }
                else
                {
                    fprintf(outf, "%d %d %d", c0, c1, c2);
                }

            }
        }
        fprintf(outf, "\n");
    }
}

int readdata(hid_t dataset, hid_t dataspace, hid_t datatype, int prj, char *data, int ntt, int npp, int nt, int np, int toff, int poff, int tmoff, int pmoff, int order, int pstride, int tstride)
{
    char buf[1001];
    hid_t memspace = -1;
    herr_t ret;
    hsize_t memdim[2], count[2], offset[2], memoffset[2], stride[2], stridecount[2], memcount[2];

    stridecount[0] = 1;
    stridecount[1] = 1;

    int ntts = ntt / tstride;
    if (ntts * tstride != ntt)ntts++;
    int nts = nt / tstride;
    if (nts * tstride != nt)nts++;
    int npps = npp / pstride;
    if (npps * pstride != npp)npps++;
    int nps = np / pstride;
    if (nps * pstride != np)nps++;

    if (order == 1)
    {
        offset[0] = toff;
        offset[1] = poff;
        count[0] = ntt;
        count[1] = npp;
        memcount[0] = ntts;
        memcount[1] = npps;
        memoffset[0] = tmoff;
        memoffset[1] = pmoff;
        memdim[0] = nts;
        memdim[1] = nps;
        stride[0] = tstride;
        stride[1] = pstride;
    }
    else
    {
        offset[1] = toff;
        offset[0] = poff;
        count[1] = ntt;
        count[0] = npp;
        memcount[1] = ntts;
        memcount[0] = npps;
        memoffset[1] = tmoff;
        memoffset[0] = pmoff;
        memdim[1] = nts;
        memdim[0] = nps;
        stride[1] = tstride;
        stride[0] = pstride;
    }
    if (tstride > 1 || pstride > 1)
    {
        ret = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, stride, memcount, stridecount);
    }
    else
    {
        ret = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);
    }
    if (ret < 0)
    {
        sprintf(buf, "Cannot select hyperslab in database\n");
        printout(buf);
        return 0;
    }
    if (tstride > 1 || pstride > 1)
    {
        memspace = H5Screate_simple(2, memdim, NULL);
        ret = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, memoffset, NULL, memcount, NULL);
    }
    else
    {
        memspace = H5Screate_simple(2, memdim, NULL);
        ret = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, memoffset, NULL, count, NULL);
    }
    if (ret < 0)
    {
        sprintf(buf, "Cannot select hyperslab in memory\n");
        printout(buf);
        return 0;
    }
    ret = H5Dread(dataset, datatype, memspace, dataspace, H5P_DEFAULT, data);
    if (ret < 0)
    {
        sprintf(buf, "Cannot read data\n");
        printout(buf);
        H5Sclose(memspace);
        return 0;
    }
    H5Sclose(memspace);

    return 1;
}

int query()
{
    char buf[1001], pvaluetype[20], tvaluetype[20], trangetype[20], prangetype[20], name[256], chrname[20], dest[10], format[10];
    char outfile[256], orientation[10], datarr[20], username[50], userpass[50];
    int iprange, itrange, ii, iformat, order, pstride, tstride, uid = -1;
    char *data = NULL;
    FILE *outf = NULL;
    long unsigned int np, *pp, *tt, nt, i;
    int prj, chr;
    hid_t project, chromosome, dataset = -1, dataspace = -1, datatype = -1;

    if (readline(username, 50, __LINE__, 1) < 0)return -1;
    if (readline(userpass, 50, __LINE__, 0) < 0)return -1;
    if (strcmp(username, "serveradmin") == 0)
    {
        if (check_pass(userpass) == 0)
        {
            sprintf(errbuf, "Invalid user name or password\n");
            gotoend();
            return 0;
        }
    }
    else
    {
        uid = check_user_pass(username, userpass);
        if (uid == -1)
        {
            sprintf(errbuf, "Invalid user name or password\n");
            gotoend();
            return 0;
        }
    }
    if (readline(pvaluetype, 20, __LINE__, 1) < 0)return -1;
    if (readline(prangetype, 20, __LINE__, 1) < 0)return -1;
    np = 0;
    if (strcmp(prangetype, "all") == 0)
    {
        np = 0;
        iprange = 2;
    }
    else if (strcmp(prangetype, "range") == 0)
    {
        np = 2;
        iprange = 1;
    }
    else
    {
        iprange = 0;
    }
    if (readline(buf, 200, __LINE__, -1) > 0)
    {
        chomp(buf);
        if (iprange == 0)
        {
            np = strtol(buf, NULL, 0);
            if (np == 0)
            {
                sprintf(errbuf, "Invalid positions number of indexes\n");
                gotoend();
                return 0;
            }
        }
    }
    else
    {
        if (iprange == 0)return -1;
    }
    if (readline(buf, 20, __LINE__, 1) < 0)return -1;
    pstride = strtol(buf, NULL, 0);
    if (pstride == 0)
    {
        sprintf(errbuf, "Invalid positions stride '%s'\n", buf);
        gotoend();
        return -1;
    }
    if (pstride != 1 && iprange == 0)
    {
        sprintf(errbuf, "Invalid positions stride '%s'\n", buf);
        gotoend();
        return -1;
    }
    if (readline(tvaluetype, 20, __LINE__, 1) < 0)return -1;
    if (readline(trangetype, 20, __LINE__, 1) < 0)return -1;
    nt = 0;
    if (strcmp(trangetype, "all") == 0)
    {
        nt = 0;
        itrange = 2;
    }
    else if (strcmp(trangetype, "range") == 0)
    {
        nt = 2;
        itrange = 1;
    }
    else
    {
        itrange = 0;
    }
    if (readline(buf, 200, __LINE__, -1) > 0)
    {
        chomp(buf);
        if (itrange == 0)
        {
            nt = strtol(buf, NULL, 0);
            if (nt == 0 && itrange == 0)
            {
                sprintf(errbuf, "Invalid taxa number of indexes\n");
                gotoend();
                return 0;
            }
        }
    }
    else
    {
        if (itrange == 0)return -1;
    }
    if (readline(buf, 20, __LINE__, 1) < 0)return -1;
    tstride = strtol(buf, NULL, 0);
    if (tstride == 0)
    {
        sprintf(errbuf, "Invalid taxa stride '%s'\n", buf);
        gotoend();
        return -1;
    }
    if (tstride != 1 && itrange == 0)
    {
        sprintf(errbuf, "Invalid taxa stride '%s'\n", buf);
        gotoend();
        return -1;
    }
    if (readline(dest, 6, __LINE__, 1) < 0)return -1;
    if (readline(format, 6, __LINE__, 1) < 0)return -1;
    iformat = 0;
    if (strcmp(format, "num") == 0)iformat = 1;
    if (readline(orientation, 9, __LINE__, 1) < 0)return -1;
    if (strcmp(orientation, "auto") != 0 && strcmp(orientation, "pf") != 0 && strcmp(orientation, "tf") != 0 && strcmp(orientation, "nodata") != 0)
    {
        sprintf(errbuf, "Invalid orientation '%s' requested\n", orientation);
        gotoend();
        return 0;
    }
    if (readline(name, 256, __LINE__, 1) < 0)return -1;
    //find project
    prj = find_project(name);
    if (prj < 0)
    {
        gotoend();
        sprintf(errbuf, "Invalid project name\n");
        return 0;
    }
    if (uid != -1)
    {
        if (accesstable[uid] == NULL)
        {
            gotoend();
            sprintf(errbuf, "Access to project denied\n");
            return 0;
        }
        sprintf(buf, "\n%s\n", name);
        if (strstr(accesstable[uid], buf) == NULL)
        {
            gotoend();
            sprintf(errbuf, "Access to project denied\n");
            return 0;
        }
    }
    if (readline(chrname, 20, __LINE__, 1) < 0)return -1;
    //find chr
    chr = find_chr(chrname, prj);
    if (chr < 0)
    {
        gotoend();
        sprintf(errbuf, "Invalid chromosome name %s\n", chrname);
        return 0;
    }

    sprintf(buf, "/mnt%d/project%d", projects[prj].mnt, projects[prj].num);
    project = H5Gopen(file, buf, H5P_DEFAULT);
    if (project < 0)
    {
        gotoend();
        printout("Cannot open project\n");
        return 0;
    }
    sprintf(buf, "chr%d", chrinfo[chr].num);
    chromosome = H5Gopen(project, buf, H5P_DEFAULT);
    if (chromosome < 0)
    {
        gotoend();
        printout("Cannot open chromosome\n");
        H5Gclose(project);
        return 0;
    }

    if (np > 0)
    {
        pp = (long unsigned int *) malloc(np * sizeof (long unsigned int));
        for (i = 0; i < np; i++)
        {
            if (readline(buf, 256, __LINE__, 1) < 0)return -1;
            if (convert_position(buf, chromosome, pp, i, pvaluetype, prangetype) < 0)
            {
                printout("cannot convert position to index\n");
                gotoend();
                free(pp);
                H5Gclose(project);
                H5Gclose(chromosome);
                return 0;
            }
        }
    }
    else
    {
        if (strcmp(prangetype, "all") == 0)
        {
            np = 2;
            pp = (long unsigned int *) malloc(np * sizeof (long unsigned int));
            pp[0] = 0;
            pp[1] = (long unsigned int) chrinfo[chr].dim0 - 1;
            //fprintf(stderr, "pos dims %lu %lu\n", pp[0], pp[1]);
        }
        else
        {
            sprintf(errbuf, "invalid range type\n");
            gotoend();
            H5Gclose(project);
            H5Gclose(chromosome);
            return 0;
        }
    }
    if (nt > 0)
    {
        tt = (long unsigned int *) malloc(nt * sizeof (long unsigned int));
        for (i = 0; i < nt; i++)
        {
            if (readline(buf, 256, __LINE__, 1) < 0)return -1;
            if (convert_taxa(buf, project, tt, i, tvaluetype) < 0)
            {
                sprintf(errbuf, "cannot convert taxa to index\n");
                gotoend();
                free(pp);
                free(tt);
                H5Gclose(project);
                H5Gclose(chromosome);
                return 0;
            }
        }
    }
    else
    {
        if (strcmp(trangetype, "all") == 0)
        {
            nt = 2;
            tt = (long unsigned int *) malloc(nt * sizeof (long unsigned int));
            tt[0] = 0;
            tt[1] = (long unsigned int) chrinfo[chr].dim1 - 1;
            //fprintf(stderr, "taxa dims %lu %lu\n", tt[0], tt[1]);
        }
        else
        {
            sprintf(errbuf, "invalid range type\n");
            free(pp);
            gotoend();
            H5Gclose(project);
            H5Gclose(chromosome);
            return 0;
        }
    }
    if (readline(buf, 256, __LINE__, 0) != -1)
    {
        printout("invalid input, expected end-of-command here\n");
        free(pp);
        free(tt);
        H5Gclose(project);
        H5Gclose(chromosome);
        return -1;
    }

    int nodata = 0;
    if (strcmp(orientation, "nodata") == 0)
    {
        nodata = 1;
    }

    if (projects[prj].order == 3)
    {
        if (strcmp(orientation, "auto") == 0)
        {
            if (iprange == 1 || iprange == 2)
            {
                if (itrange == 1 || itrange == 2)
                {
                    if (pp[1] - pp[0] >= tt[1] - tt[0])
                    {
                        order = 1;
                    }
                    else
                    {
                        order = 2;
                    }
                }
                else
                {
                    order = 1;
                }
            }
            else
            {
                if (itrange == 1 || itrange == 2)
                {
                    order = 2;
                }
                else
                {
                    order = 1;
                }
            }
        }
        else if (strcmp(orientation, "pf") == 0)
        {
            order = 1;
        }
        else
        {
            order = 2;
        }
    }
    else
    {
        order = projects[prj].order;
    }
    if (order == 1)
    {
        sprintf(datarr, "%s_pf", DATAARRAY);
    }
    else
    {
        sprintf(datarr, "%s_tf", DATAARRAY);
    }

    outf = outs;
    outfile[0] = '\0';
    if (strcmp(dest, "file") == 0)
    {
        unsigned int iseed = (unsigned int) time(NULL);
        srand(iseed);

        sprintf(outfile, "%d_%s_%s_%u%u", projects[prj].mnt, projects[prj].name, chrinfo[chr].name, iseed, rand());
        sprintf(prtbuf, "%s/%s", tmpdir, outfile);
        if ((outf = fopen(prtbuf, "w")) == NULL)
        {
            sprintf(buf, "Cannot open output file %s\n", prtbuf);
            printout(buf);
            free(pp);
            free(tt);
            H5Gclose(project);
            H5Gclose(chromosome);
            return 0;
        }
        //fprintf(outs, "%s\n", outfile);
    }

    //output data start
    if (outfile[0] != '\0')
    {
        sprintf(prtbuf, "%s\n", outfile);
        printout(prtbuf);
    }
    fprintf(outf, "encoding = %s\n", projects[prj].encoding);
    fprintf(outf, "data element size %d bytes\n", projects[prj].enc);
    if (order == 1)
    {
        fprintf(outf, "orientation = pf\n");
    }
    else
    {
        fprintf(outf, "orientation = tf\n");
    }

    if (nodata == 0)
    {
        dataset = H5Dopen(chromosome, datarr, H5P_DEFAULT);
        if (dataset < 0)
        {
            sprintf(buf, "Cannot open dataset '%s'\n", datarr);
            printout(buf);
            free(pp);
            free(tt);
            H5Gclose(project);
            H5Gclose(chromosome);
            return 0;
        }
        dataspace = H5Dget_space(dataset);
        if (dataspace < 0)
        {
            sprintf(buf, "Cannot open dataspace in dataset '%s'\n", datarr);
            printout(buf);
            free(pp);
            free(tt);
            H5Gclose(project);
            H5Gclose(chromosome);
            H5Dclose(dataset);
            return 0;
        }
        datatype = H5Dget_type(dataset);
        if (datatype < 0)
        {
            sprintf(buf, "Cannot get datatype of dataset '%s'\n", datarr);
            free(pp);
            free(tt);
            H5Gclose(project);
            H5Gclose(chromosome);
            H5Dclose(dataset);
            H5Sclose(dataspace);
            printout(buf);
            return 0;
        }
    }

    if (iprange == 1 || iprange == 2)
    {
        if (itrange == 1 || itrange == 2)
        {
            //range-range
            int npp = pp[1] - pp[0] + 1;
            int ntt = tt[1] - tt[0] + 1;
            int npps = npp / pstride;
            int ntts = ntt / tstride;
            long unsigned int memory = (long unsigned int)npps * (long unsigned int)ntts * (long unsigned int)projects[prj].enc * sizeof (char);
            if (npps * pstride != npp)npps++;
            if (ntts * tstride != ntt)ntts++;
            fprintf(outf, "positions size = %d\n", npps);
            fprintf(outf, "taxa size = %d\n", ntts);

            if (nodata == 0)
            {
                data = malloc(memory);
                if (data == NULL)
                {
                    sprintf(buf, "Cannot allocate %lu bytes of memory [range-range]\n", memory);
                    printout(buf);
                    free(pp);
                    free(tt);
                    return 0;
                }

                if (readdata(dataset, dataspace, datatype, prj, data, ntt, npp, ntt, npp, tt[0], pp[0], 0, 0, order, pstride, tstride) <= 0)
                {
                    free(pp);
                    free(tt);
                    free(data);
                    H5Gclose(project);
                    H5Gclose(chromosome);
                    H5Dclose(dataset);
                    H5Tclose(datatype);
                    H5Sclose(dataspace);
                    return 0;
                }
                printdata(prj, data, iformat, ntt, npp, outf, pstride, tstride, order);
            }
            printtable("positions", chrname, project, pp[0], pp[1], prj, pstride, outf);
            printtable("alleles", chrname, project, pp[0], pp[1], prj, pstride, outf);
            printtable("taxa", chrname, project, tt[0], tt[1], prj, tstride, outf);
            printtable("markers", chrname, project, pp[0], pp[1], prj, pstride, outf);
        }
        else
        {
            //range-list
            int npp = pp[1] - pp[0] + 1;
            int npps = npp / pstride;
            if (npps * pstride != npp)npps++;
            int ntt = 1;
            fprintf(outf, "positions size = %i\n", npps);
            fprintf(outf, "taxa size = %li\n", nt);
            
            long unsigned int memory = (long unsigned int)npps * (long unsigned int)nt * (long unsigned int)projects[prj].enc * sizeof (char);

            if (nodata == 0)
            {
                data = malloc(memory);
                if (data == NULL)
                {
                    sprintf(buf, "Cannot allocate %lu bytes of memory [range-list]\n", memory);
                    printout(buf);
                    H5Gclose(project);
                    H5Gclose(chromosome);
                    H5Dclose(dataset);
                    H5Tclose(datatype);
                    H5Sclose(dataspace);
                    free(pp);
                    free(tt);
                    return 0;
                }

                for (ii = 0; ii < nt; ii++)
                {
                    if (readdata(dataset, dataspace, datatype, prj, data, ntt, npp, nt, npp, tt[ii], pp[0], ii, 0, order, pstride, tstride) <= 0)
                    {
                        free(pp);
                        free(tt);
                        free(data);
                        H5Gclose(project);
                        H5Gclose(chromosome);
                        H5Dclose(dataset);
                        H5Tclose(datatype);
                        H5Sclose(dataspace);
                        return 0;
                    }
                }
                printdata(prj, data, iformat, nt, npp, outf, pstride, tstride, order);
            }
            printtable("positions", chrname, project, pp[0], pp[1], prj, pstride, outf);
            printtable("alleles", chrname, project, pp[0], pp[1], prj, pstride, outf);
            for (ii = 0; ii < nt; ii++)
                printtable("taxa", chrname, project, tt[ii], tt[ii], prj, 1, outf);
            printtable("markers", chrname, project, pp[0], pp[1], prj, pstride, outf);
        }
    }
    else
    {
        if (itrange == 1 || itrange == 2)
        {
            //list-range
            int npp = 1;
            int ntt = tt[1] - tt[0] + 1;
            int ntts = ntt / tstride;
            if (ntts * tstride != ntt)ntts++;
            fprintf(outf, "positions size = %li\n", np);
            fprintf(outf, "taxa size = %i\n", ntts);
            
            long unsigned int memory = (long unsigned int)np * (long unsigned int)ntts * (long unsigned int)projects[prj].enc * sizeof (char);

            if (nodata == 0)
            {
                data = malloc(memory);
                if (data == NULL)
                {
                    sprintf(buf, "Cannot allocate %lu bytes of memory [list-range]\n", memory);
                    printout(buf);
                    H5Gclose(project);
                    H5Gclose(chromosome);
                    H5Dclose(dataset);
                    H5Tclose(datatype);
                    H5Sclose(dataspace);
                    free(pp);
                    free(tt);
                    return 0;
                }
                for (ii = 0; ii < np; ii++)
                {
                    if (readdata(dataset, dataspace, datatype, prj, data, ntt, npp, ntt, np, tt[0], pp[ii], 0, ii, order, pstride, tstride) <= 0)
                    {
                        free(pp);
                        free(tt);
                        free(data);
                        H5Gclose(project);
                        H5Gclose(chromosome);
                        H5Dclose(dataset);
                        H5Tclose(datatype);
                        H5Sclose(dataspace);
                        return 0;
                    }
                }
                printdata(prj, data, iformat, ntt, np, outf, pstride, tstride, order);
            }
            for (ii = 0; ii < np; ii++)
                printtable("positions", chrname, project, pp[ii], pp[ii], prj, 1, outf);
            for (ii = 0; ii < np; ii++)
                printtable("alleles", chrname, project, pp[ii], pp[ii], prj, 1, outf);
            printtable("taxa", chrname, project, tt[0], tt[1], prj, tstride, outf);
            for (ii = 0; ii < np; ii++)
                printtable("markers", chrname, project, pp[ii], pp[ii], prj, 1, outf);
        }
        else
        {
            //list-list
            if (np != nt)
            {
                sprintf(buf, "[list-list]: taxa and positions numbers do not match\n");
                printout(buf);
                H5Gclose(project);
                H5Gclose(chromosome);
                H5Dclose(dataset);
                H5Tclose(datatype);
                H5Sclose(dataspace);
                free(pp);
                free(tt);
                return 0;
            }
            int npp = 1;
            int ntt = 1;
            fprintf(outf, "positions size = %li\n", np);
            fprintf(outf, "taxa size = %li\n", nt);

            if (nodata == 0)
            {
                data = malloc((long unsigned int)projects[prj].enc * sizeof (char));
                if (data == NULL)
                {
                    sprintf(buf, "Cannot allocate %lu bytes of memory [list-list]\n", (long unsigned int)projects[prj].enc * sizeof (char));
                    printout(buf);
                    H5Gclose(project);
                    H5Gclose(chromosome);
                    H5Dclose(dataset);
                    H5Tclose(datatype);
                    H5Sclose(dataspace);
                    free(pp);
                    free(tt);
                    return 0;
                }
                int kk = 1;
                for (ii = 0; ii < nt; ii++)
                {
                    if (readdata(dataset, dataspace, datatype, prj, data, ntt, npp, ntt, npp, tt[ii], pp[ii], 0, 0, order, pstride, tstride) <= 0)
                    {
                        H5Gclose(project);
                        H5Gclose(chromosome);
                        H5Dclose(dataset);
                        H5Tclose(datatype);
                        H5Sclose(dataspace);
                        free(pp);
                        free(tt);
                        free(data);
                        return 0;
                    }
                    printdata(prj, data, iformat, ntt, npp, outf, pstride, tstride, order);
                    kk++;
                }
            }
            for (ii = 0; ii < np; ii++)
                printtable("positions", chrname, project, pp[ii], pp[ii], prj, 1, outf);
            for (ii = 0; ii < np; ii++)
                printtable("alleles", chrname, project, pp[ii], pp[ii], prj, 1, outf);
            for (ii = 0; ii < nt; ii++)
                printtable("taxa", chrname, project, tt[ii], tt[ii], prj, 1, outf);
            for (ii = 0; ii < np; ii++)
                printtable("markers", chrname, project, pp[ii], pp[ii], prj, 1, outf);
        }
    }
    //free memory
    free(pp);
    free(tt);
    free(data);
    H5Gclose(project);
    H5Gclose(chromosome);
    if (nodata == 0)
    {
        H5Dclose(dataset);
        H5Tclose(datatype);
        H5Sclose(dataspace);
    }


    if (strcmp(dest, "file") == 0)fclose(outf);
    return 1;
}

int quit()
{
    char pass[100];
    char buf1[257];

    if (readline(pass, 100, __LINE__, 0) < 0)return -1;
    if (readline(buf1, 256, __LINE__, 0) != -1)
    {
        printout("invalid input, expected end-of-command here\n");
        return -1;
    }
    if (check_pass(pass) == 0)
    {
        printout("invalid password\n");
        return 0;
    }
    shutdown();
    return 1;
}

int cache()
{
    char buf1[257];
    int i;

    if (readline(buf1, 256, __LINE__, 0) != -1)
    {
        printout("invalid input, expected end-of-command here\n");
        return -1;
    }

    fprintf(outs, "%d files\n", nfiles);
    for (i = 0; i < nfiles; i++)
    {
        fprintf(outs, "/mnt%d\t%s\t\n", files[i].mnt, files[i].name);
    }
    fprintf(outs, "%d projects\n", nprojects);
    for (i = 0; i < nprojects; i++)
    {
        char ordr[10];
        if (projects[i].order == 1)strcpy(ordr, "pf");
        if (projects[i].order == 2)strcpy(ordr, "tf");
        if (projects[i].order == 3)strcpy(ordr, "both");
        fprintf(outs, "/mnt%d/project%d\t%s\tenc_size=%dB\torientation=%s\n", projects[i].mnt, projects[i].num, projects[i].name, projects[i].enc, ordr);
    }
    fprintf(outs, "%d chromosomes\n", nchrinfo);
    for (i = 0; i < nchrinfo; i++)
    {
        fprintf(outs, "/mnt%d/project%d/chr%d\t%s\tproject_indx=%d\tpos_dim=%lu\ttaxa_dim=%lu\n", chrinfo[i].mnt, projects[chrinfo[i].prj].num, chrinfo[i].num, chrinfo[i].name, chrinfo[i].prj, (long unsigned int) chrinfo[i].dim0, (long unsigned int) chrinfo[i].dim1);
    }
    fflush(outs);
    return 1;
}

int writestrtable(), writeinttable();

int verifyorcreate(hid_t obj, char *name, hsize_t *mcount)
{
    hid_t dataset, dataspace;
    hsize_t dims[2];
    herr_t ret;
    hsize_t mastercount;
    char namei[256], namev[256];


    sprintf(namei, "%s_si", name);
    sprintf(namev, "%s_sv", name);
    dataset = H5Dopen(obj, name, H5P_DEFAULT);
    if (dataset < 0)
    {
        sprintf(prtbuf, "Cannot open dataset '%s'\n", name);
        printout(prtbuf);
        return -1;
    }
    dataspace = H5Dget_space(dataset);
    ret = H5Sget_simple_extent_dims(dataspace, dims, NULL);
    mastercount = dims[0];
    *mcount = dims[0];
    H5Dclose(dataset);
    H5Sclose(dataspace);

    if (H5Lexists(obj, namev, H5P_DEFAULT))
    {
        ret = H5Ldelete(obj, namev, H5P_DEFAULT);
        H5Fflush(file, H5F_SCOPE_LOCAL);
        if (ret < 0)
        {
            sprintf(prtbuf, "Cannot delete dataset '%s'\n", namev);
            printout(prtbuf);
            return -1;
        }
    }

    if (H5Lexists(obj, namei, H5P_DEFAULT))
    {
        ret = H5Ldelete(obj, namei, H5P_DEFAULT);
        H5Fflush(file, H5F_SCOPE_LOCAL);
        if (ret < 0)
        {
            sprintf(prtbuf, "Cannot delete dataset '%s'\n", namei);
            printout(prtbuf);
            return -1;
        }
    }
    return 1;
}

int read_access_table()
{
    char buf[1001];
    hid_t dataspace, dataset, datatype, xfer_pid;
    herr_t ret;
    hsize_t dims[2];
    char name[10];

    strcpy(name, "/access");

    dataset = H5Dopen(file, name, H5P_DEFAULT);
    if (dataset < 0)
    {
        sprintf(buf, "Cannot open dataset '%s'\n", name);
        printout(buf);
        return -1;
    }
    dataspace = H5Dget_space(dataset);
    if (dataspace < 0)
    {
        sprintf(buf, "Cannot open dataspace in dataset '%s'\n", name);
        printout(buf);
        H5Dclose(dataset);
        return 0;
    }
    datatype = H5Dget_type(dataset);
    if (datatype < 0)
    {
        sprintf(buf, "Cannot get datatype of dataset '%s'\n", name);
        H5Dclose(dataset);
        H5Dclose(dataspace);
        printout(buf);
        return -1;
    }

    dims[0] = nusers;

    accesstable = (char **) malloc(dims[0] * sizeof (char *));
    if (accesstable == NULL)
    {
        sprintf(buf, "Cannot get memory for dataset '%s'\n", name);
        printout(buf);
        H5Dclose(dataset);
        H5Dclose(dataspace);
        H5Tclose(datatype);
        return -1;
    }

    //fprintf(stderr, "mem allocated for strtable (%llu)\n", size);

    xfer_pid = H5Pcreate(H5P_DATASET_XFER);
    if (xfer_pid < 0)
    {
        sprintf(buf, "Cannot get xfer_pid of dataset '%s'\n", name);
        printout(buf);
        H5Dclose(dataset);
        H5Tclose(datatype);
        H5Dclose(dataspace);
        free(accesstable);
        return -1;
    }

    //fprintf(stderr, "ready to read\n");
    ret = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, xfer_pid, accesstable);
    if (ret < 0)
    {
        sprintf(buf, "Cannot read data from %s\n", name);
        printout(buf);
        H5Dclose(dataset);
        H5Sclose(dataspace);
        H5Tclose(datatype);
        H5Pclose(xfer_pid);
        free(accesstable);
        return -1;
    }
    //fprintf(stderr, "read done\n");

    H5Tclose(datatype);
    H5Dclose(dataset);
    H5Sclose(dataspace);
    H5Pclose(xfer_pid);

    return 1;
}

int readstrtable(arrsort **sortstrtableout, char ***strtableout, hid_t obj, char *name, hsize_t size)
{
    char buf[1001];
    hid_t dataspace, dataset, datatype, xfer_pid;
    herr_t ret;
    hsize_t dims[2], i;
    char **strtable;
    arrsort *sortstrtable;

    dataset = H5Dopen(obj, name, H5P_DEFAULT);
    if (dataset < 0)
    {
        sprintf(buf, "Cannot open dataset '%s'\n", name);
        printout(buf);
        return -1;
    }
    dataspace = H5Dget_space(dataset);
    if (dataspace < 0)
    {
        sprintf(buf, "Cannot open dataspace in dataset '%s'\n", name);
        printout(buf);
        H5Dclose(dataset);
        return 0;
    }
    datatype = H5Dget_type(dataset);
    if (datatype < 0)
    {
        sprintf(buf, "Cannot get datatype of dataset '%s'\n", name);
        H5Dclose(dataset);
        H5Dclose(dataspace);
        printout(buf);
        return -1;
    }

    dims[0] = size;

    strtable = (char **) malloc(dims[0] * sizeof (char *));
    if (strtable == NULL)
    {
        sprintf(buf, "Cannot get memory for dataset '%s'\n", name);
        printout(buf);
        H5Dclose(dataset);
        H5Dclose(dataspace);
        H5Tclose(datatype);
        return -1;
    }

    //fprintf(stderr, "mem allocated for strtable (%llu)\n", size);

    xfer_pid = H5Pcreate(H5P_DATASET_XFER);
    if (xfer_pid < 0)
    {
        sprintf(buf, "Cannot get xfer_pid of dataset '%s'\n", name);
        printout(buf);
        H5Dclose(dataset);
        H5Tclose(datatype);
        H5Dclose(dataspace);
        free(strtable);
        return -1;
    }

    //fprintf(stderr, "ready to read\n");
    ret = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, xfer_pid, strtable);
    if (ret < 0)
    {
        sprintf(buf, "Cannot read data from %s\n", name);
        printout(buf);
        H5Dclose(dataset);
        H5Sclose(dataspace);
        H5Tclose(datatype);
        H5Pclose(xfer_pid);
        free(strtable);
        return -1;
    }
    //fprintf(stderr, "read done\n");

    H5Tclose(datatype);
    H5Dclose(dataset);
    H5Sclose(dataspace);
    H5Pclose(xfer_pid);

    sortstrtable = (arrsort *) malloc(dims[0] * sizeof (arrsort));
    if (sortstrtable == NULL)
    {
        sprintf(buf, "Cannot get memory (2) for dataset '%s'\n", name);
        printout(buf);
        free(strtable);
        return -1;
    }

    //fprintf(stderr, "mem allocated for sortstrtable (%llu)\n", size);

    for (i = 0; i < dims[0]; i++)
    {
        //fprintf(stderr, "copy %llu int\n", i);
        sortstrtable[i].index = (long unsigned int) i;
        //fprintf(stderr, "copy %llu str %s\n", i, strtable[i]);
        sortstrtable[i].string = strtable[i];
    }

    *sortstrtableout = sortstrtable;
    *strtableout = strtable;
    return 1;
}

void freestrtable1(char **strtable, unsigned int size, unsigned int ex)
{
    unsigned int i;

    for (i = 0; i < size; i++)
    {
        if (i != ex)free(strtable[i]);
    }
    free(strtable);
}

void freestrtable(char **strtable, unsigned int size)
{
    unsigned int i;

    for (i = 0; i < size; i++)
    {
        free(strtable[i]);
    }
    free(strtable);
}

int arrsort_cmp(const void *a, const void *b)
{
    const arrsort *ia = (const arrsort *) a; // casting pointer types 
    const arrsort *ib = (const arrsort *) b;
    return strcmp(ia->string, ib->string);
    /* returns negative if b > a 
    and positive if a > b */
}

int dosortstrtable(char **strtable, arrsort *sortstrtable, int **inttableout, unsigned int size)
{
    unsigned int i;
    int *inttable;

    //fprintf(stderr, "Ready to do qsort\n");

    qsort(sortstrtable, size, sizeof (arrsort), arrsort_cmp);

    //fprintf(stderr, "qsort done\n");

    //create new strtable and inttable

    inttable = (int *) malloc(size * sizeof (int));
    if (inttable == NULL)
    {
        sprintf(prtbuf, "Cannot get memory (4) for int table\n");
        printout(prtbuf);
        freestrtable(strtable, size);
        free(sortstrtable);
        return -1;
    }

    //fprintf(stderr, "got memory for inttable\n");
    for (i = 0; i < size; i++)
    {
        strtable[i] = sortstrtable[i].string;
        inttable[i] = (unsigned int) sortstrtable[i].index;
    }

    *inttableout = inttable;
    return 1;
}

int write_access_table() // (char **strtable, hid_t obj, char *nameorig, hsize_t size)
{
    char buf[1001], name[266];
    hid_t dataspace, dataset, datatype, xfer_pid, cparms;
    herr_t ret;
    hsize_t dims[2];

    sprintf(name, "/access");

    dims[0] = nusers;
    if (nusers == 0)dims[0] = 1;

    //fprintf(stderr, "write %s size %llu\n", name, size);

    dataspace = H5Screate_simple(1, dims, NULL);
    if (dataspace < 0)
    {
        sprintf(buf, "Cannot create dataspace in dataset '%s'\n", name);
        printout(buf);
        return 0;
    }
    datatype = H5Tcopy(H5T_C_S1);
    if (datatype < 0)
    {
        sprintf(buf, "Cannot get datatype of dataset '%s'\n", name);
        H5Dclose(dataspace);
        printout(buf);
        return -1;
    }
    ret = H5Tset_size(datatype, H5T_VARIABLE);
    cparms = H5Pcreate(H5P_DATASET_CREATE);
    dataset = H5Dcreate(file, name, datatype, dataspace, H5P_DEFAULT, cparms, H5P_DEFAULT);
    if (dataset < 0)
    {
        sprintf(buf, "Cannot create dataset '%s'\n", name);
        H5Tclose(datatype);
        H5Dclose(dataspace);
        printout(buf);
        return -1;
    }

    xfer_pid = H5Pcreate(H5P_DATASET_XFER);
    if (xfer_pid < 0)
    {
        sprintf(buf, "Cannot get xfer_pid of dataset '%s'\n", name);
        H5Dclose(dataset);
        H5Tclose(datatype);
        H5Dclose(dataspace);
        printout(buf);
        return -1;
    }

    ret = H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, xfer_pid, accesstable);
    if (ret < 0)
    {
        sprintf(buf, "Cannot write data to %s\n", name);
        printout(buf);
        H5Dclose(dataset);
        H5Sclose(dataspace);
        H5Tclose(datatype);
        return -1;
    }
    H5Tclose(datatype);
    H5Dclose(dataset);
    H5Sclose(dataspace);
    H5Fflush(file, H5F_SCOPE_LOCAL);
    return 1;
}

int writestrtable(char **strtable, hid_t obj, char *nameorig, hsize_t size)
{
    char buf[1001], name[266];
    hid_t dataspace, dataset, datatype, xfer_pid, cparms;
    herr_t ret;
    hsize_t dims[2];

    sprintf(name, "%s_sv", nameorig);

    dims[0] = size;

    //fprintf(stderr, "write %s size %llu\n", name, size);

    dataspace = H5Screate_simple(1, dims, NULL);
    if (dataspace < 0)
    {
        sprintf(buf, "Cannot create dataspace in dataset '%s'\n", name);
        printout(buf);
        return 0;
    }
    datatype = H5Tcopy(H5T_C_S1);
    if (datatype < 0)
    {
        sprintf(buf, "Cannot get datatype of dataset '%s'\n", name);
        H5Dclose(dataspace);
        printout(buf);
        return -1;
    }
    ret = H5Tset_size(datatype, H5T_VARIABLE);
    cparms = H5Pcreate(H5P_DATASET_CREATE);
    dataset = H5Dcreate(obj, name, datatype, dataspace, H5P_DEFAULT, cparms, H5P_DEFAULT);
    if (dataset < 0)
    {
        sprintf(buf, "Cannot create dataset '%s'\n", name);
        H5Tclose(datatype);
        H5Dclose(dataspace);
        printout(buf);
        return -1;
    }

    xfer_pid = H5Pcreate(H5P_DATASET_XFER);
    if (xfer_pid < 0)
    {
        sprintf(buf, "Cannot get xfer_pid of dataset '%s'\n", name);
        H5Dclose(dataset);
        H5Tclose(datatype);
        H5Dclose(dataspace);
        printout(buf);
        return -1;
    }

    ret = H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, xfer_pid, strtable);
    if (ret < 0)
    {
        sprintf(buf, "Cannot write data to %s\n", name);
        printout(buf);
        H5Dclose(dataset);
        H5Sclose(dataspace);
        H5Tclose(datatype);
        return -1;
    }
    H5Tclose(datatype);
    H5Dclose(dataset);
    H5Sclose(dataspace);
    H5Fflush(file, H5F_SCOPE_LOCAL);
    return 1;
}

int writeinttable(int *inttable, hid_t obj, char *nameorig, hsize_t size)
{
    char buf[1001], name[266];
    hid_t dataspace, dataset, datatype, xfer_pid, cparms;
    herr_t ret;
    hsize_t dims[2];

    sprintf(name, "%s_si", nameorig);

    dims[0] = size;
    dataspace = H5Screate_simple(1, dims, NULL);
    if (dataspace < 0)
    {
        sprintf(buf, "Cannot open dataspace in dataset '%s'\n", name);
        printout(buf);
        return 0;
    }
    datatype = H5Tcopy(H5T_NATIVE_INT);
    cparms = H5Pcreate(H5P_DATASET_CREATE);
    dataset = H5Dcreate(obj, name, datatype, dataspace, H5P_DEFAULT, cparms, H5P_DEFAULT);
    if (dataset < 0)
    {
        sprintf(buf, "Cannot open dataset '%s'\n", name);
        H5Dclose(dataspace);
        H5Tclose(datatype);
        printout(buf);
        return -1;
    }

    xfer_pid = H5Pcreate(H5P_DATASET_XFER);
    if (xfer_pid < 0)
    {
        sprintf(buf, "Cannot get xfer_pid of dataset '%s'\n", name);
        H5Dclose(dataset);
        H5Tclose(datatype);
        H5Dclose(dataspace);
        printout(buf);
        return -1;
    }

    ret = H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, xfer_pid, inttable);
    if (ret < 0)
    {
        sprintf(buf, "Cannot write data to %s\n", name);
        printout(buf);
        H5Dclose(dataset);
        H5Sclose(dataspace);
        H5Tclose(datatype);
        return -1;
    }
    H5Tclose(datatype);
    H5Dclose(dataset);
    H5Sclose(dataspace);
    H5Fflush(file, H5F_SCOPE_LOCAL);
    return 1;
}

int writetables(char **strtable, int *inttable, hid_t obj, char *nameorig, hsize_t size)
{
    if (writestrtable(strtable, obj, nameorig, size) == -1)return -1;
    return writeinttable(inttable, obj, nameorig, size);
}

int sort()
{
    char buf[1001], name[256], pass[100];
    hid_t project = -1, chromosome = -1, locfile;
    int i, nn;
    char buf1[257];
    char **strtable;
    int *inttable;
    arrsort *sortstrtable;
    hsize_t mcount;
    herr_t ret;


    if (readline(pass, 100, __LINE__, 0) < 0)return -1;
    if (readline(name, 255, __LINE__, 1) < 0)return -1;
    if (readline(buf1, 256, __LINE__, 0) != -1)
    {
        printout("invalid input, expected end-of-command here\n");
        return -1;
    }
    if (check_pass(pass) == 0)
    {
        printout("invalid password\n");
        return 0;
    }

    nn = find_file(name);
    if (nn < 0)
    {
        printout("invalid file name\n");
        return 0;
    }
    //umnount and close the file
    sprintf(buf, "/mnt%d", nn);
    ret = H5Funmount(file, buf);
    if (ret < 0)
    {
        printout("Cannot unmount file\n");
        return 0;
    }
    ret = H5Fclose(filesh[nn]);
    if (ret < 0)
    {
        printout("Cannot close file\n");
        if (H5Fmount(file, buf, filesh[nn], H5P_DEFAULT) < 0)
        {
            sprintf(prtbuf, "Cannot re-mount file %s on group %s", name, buf);
            error(prtbuf, __LINE__);
        }
        return 0;
    }
    sprintf(buf, "%s%s", hdf5dir, files[nn].name);
    locfile = H5Fopen(buf, H5F_ACC_RDWR, H5P_DEFAULT);
    if (locfile < 0)
    {
        printout("Cannot reopen file for writing\n");
        sprintf(buf, "%s%s", hdf5dir, files[nn].name);
        filesh[nn] = H5Fopen(buf, H5P_DEFAULT, H5P_DEFAULT);
        if (filesh[nn] < 0)
        {
            sprintf(prtbuf, "Cannot reopen file\n");
            error(prtbuf, __LINE__);
        }
        //remount
        sprintf(buf, "/mnt%d", files[nn].mnt);
        if (H5Fmount(file, buf, filesh[nn], H5P_DEFAULT) < 0)
        {
            sprintf(prtbuf, "Cannot mount file %s on group %s", name, buf);
            error(prtbuf, __LINE__);
        }
        return 0;
    }

    for (i = 0; i < nprojects; i++)
    {
        if (projects[i].mnt == nn)
        {
            //open project
            sprintf(buf, "/project%d", projects[i].num);
            project = H5Gopen(locfile, buf, H5P_DEFAULT);
            if (project < 0)
            {
                sprintf(prtbuf, "Cannot open %s in file %d\n", buf, nn);
                printout(prtbuf);
                return 0;
            }
            //sprintf(prtbuf, "project %d verifyorcreate: taxa\n", i);
            //printout(prtbuf);
            if (verifyorcreate(project, "taxa", &mcount) == -1)
            {
                H5Gclose(project);
                return 0;
            }
            //read in taxa table'
            //sprintf(prtbuf, "project %d readstrtable: taxa\n", i);
            //printout(prtbuf);
            if (readstrtable(&sortstrtable, &strtable, project, "taxa", mcount) == -1)
            {
                H5Gclose(project);
                return 0;
            }
            //sort it
            //sprintf(prtbuf, "project %d dosortstrtable: taxa\n", i);
            //printout(prtbuf);
            if (dosortstrtable(strtable, sortstrtable, &inttable, mcount) == -1)
            {
                H5Gclose(project);
                free(sortstrtable);
                freestrtable(strtable, mcount);
                return 0;
            }
            //write back sorted values to file
            //sprintf(prtbuf, "project %d writetables: taxa\n", i);
            //printout(prtbuf);
            if (writetables(strtable, inttable, project, "taxa", mcount) == -1)
            {
                H5Gclose(project);
                free(inttable);
                free(sortstrtable);
                freestrtable(strtable, mcount);
                return 0;
            }
            H5Gclose(project);
            free(inttable);
            free(sortstrtable);
            freestrtable(strtable, mcount);
        }
    }

    for (i = 0; i < nchrinfo; i++)
    {
        if (chrinfo[i].mnt == nn)
        {
            //open chromosome
            sprintf(buf, "/project%d/chr%d", projects[chrinfo[i].prj].num, chrinfo[i].num);
            chromosome = H5Gopen(locfile, buf, H5P_DEFAULT);
            if (project < 0)
            {
                sprintf(prtbuf, "Cannot open %s in file %d\n", buf, nn);
                return 0;
            }
            if (H5Lexists(chromosome, "markers", H5P_DEFAULT))
            {
                if (verifyorcreate(chromosome, "markers", &mcount) == -1)
                {
                    H5Gclose(chromosome);
                    return -1;
                }
                if (readstrtable(&sortstrtable, &strtable, chromosome, "markers", mcount) == -1)
                {
                    H5Gclose(chromosome);
                    return -1;
                }
                //sort it
                if (dosortstrtable(strtable, sortstrtable, &inttable, mcount) == -1)
                {
                    H5Gclose(chromosome);
                    free(sortstrtable);
                    freestrtable(strtable, mcount);
                    return -1;
                }
                //write back sorted values to file
                if (writetables(strtable, inttable, chromosome, "markers", mcount) == -1)
                {
                    H5Gclose(chromosome);
                    free(inttable);
                    free(sortstrtable);
                    freestrtable(strtable, mcount);
                    return -1;
                }
                free(inttable);
                free(sortstrtable);
                freestrtable(strtable, mcount);
            }
            H5Gclose(chromosome);
        }
    }

    //close local file
    ret = H5Fclose(locfile);
    if (ret < 0)
    {
        printout("Cannot close local file\n");
        return 0;
    }
    sprintf(buf, "%s%s", hdf5dir, files[nn].name);
    filesh[nn] = H5Fopen(buf, H5P_DEFAULT, H5P_DEFAULT);
    if (filesh[nn] < 0)
    {
        sprintf(prtbuf, "Cannot reopen file\n");
        error(prtbuf, __LINE__);
    }
    //remount
    sprintf(buf, "/mnt%d", files[nn].mnt);
    if (H5Fmount(file, buf, filesh[nn], H5P_DEFAULT) < 0)
    {
        sprintf(prtbuf, "Cannot mount file %s on group %s", name, buf);
        error(prtbuf, __LINE__);
    }
    return 1;
}

int table()
{
    char buf[1001], name[256], projectname[100], chrname[100], username[50], userpass[50], postype[100], indx0[100], indx1[100];
    hid_t project, chromosome;
    unsigned int n0, n1;
    char buf1[257];
    int uid = -1, prj, i, chr;

    if (readline(username, 50, __LINE__, 1) < 0)return -1;
    if (readline(userpass, 50, __LINE__, 0) < 0)return -1;
    if (strcmp(username, "serveradmin") == 0)
    {
        if (check_pass(userpass) == 0)
        {
            sprintf(errbuf, "Invalid user name or password\n");
            gotoend();
            return 0;
        }
    }
    else
    {
        uid = check_user_pass(username, userpass);
        if (uid == -1)
        {
            sprintf(errbuf, "Invalid user name or password\n");
            gotoend();
            return 0;
        }
    }
    if (readline(name, 100, __LINE__, 1) < 0)return -1;
    if (readline(projectname, 255, __LINE__, 1) < 0)return -1;
    if (readline(chrname, 100, __LINE__, -1) < 0)
    {
        if (strcmp(name, "taxa") != 0)
        {
            return -1;
        }
    }
    chomp(chrname);
    strcpy(postype, "index");
    if (readline(postype, 100, __LINE__, -1) < 0)
    {
        if (strcmp(name, "positions") != 0)
        {
            return -1;
        }
    }
    if (readline(indx0, 100, __LINE__, 1) < 0)return -1;
    n0 = strtol(indx0, NULL, 0);
    if (n0 == 0 && strcmp(indx0, "0") != 0)
    {
        gotoend();
        printout("Invalid starting index\n");
        return -1;
    }
    if (readline(indx1, 100, __LINE__, 1) < 0)return -1;
    n1 = strtol(indx1, NULL, 0);
    if ((n1 == 0 && strcmp(indx1, "0") != 0) || n1 < n0)
    {
        gotoend();
        printout("Invalid ending index\n");
        return -1;
    }

    if (readline(buf1, 256, __LINE__, 0) != -1)
    {
        printout("invalid input, expected end-of-command here\n");
        return -1;
    }

    //open project
    prj = find_project(projectname);
    if (prj < 0)
    {
        printout("Invalid project name\n");
        return 0;
    }

    if (uid != -1)
    {
        if (accesstable[uid] == NULL)
        {
            printout("Access to project denied\n");
            return 0;
        }
        sprintf(buf1, "\n%s\n", projectname);
        if (strstr(accesstable[uid], buf1) == NULL)
        {
            printout("Access to project denied\n");
            return 0;
        }
    }

    sprintf(buf, "/mnt%d/project%d", projects[prj].mnt, projects[prj].num);
    project = H5Gopen(file, buf, H5P_DEFAULT);
    if (project < 0)
    {
        printout("Cannot open project\n");
        return 0;
    }

    if (strcmp(name, "positions") == 0 && strcmp(postype, "value") == 0)
    {
        chr = find_chr(chrname, prj);
        if (chr < 0)
        {
            printout("Invalid chromosome name\n");
            return 0;
        }
        sprintf(buf, "chr%d", chrinfo[chr].num);
        chromosome = H5Gopen(project, buf, H5P_DEFAULT);
        if (chromosome < 0)
        {
            printout("Cannot open chromosome\n");
            H5Gclose(project);
            return 0;
        }
        long unsigned int pp[2];
        for (i = 0; i < 2; i++)
        {
            if (convert_position(buf, chromosome, pp, i, "positions", "range") < 0)
            {
                printout("cannot convert position to index\n");
                H5Gclose(project);
                H5Gclose(chromosome);
                return 0;
            }
        }
        H5Gclose(chromosome);
        n0 = (unsigned int) pp[0];
        n1 = (unsigned int) pp[1];
    }

    fprintf(outs, "%d\n", n0);
    fprintf(outs, "%d\n", n1);
    int retcode = printtable(name, chrname, project, n0, n1, prj, 1, outs);

    if (retcode == -1)
    {
        H5Gclose(project);
        sprintf(prtbuf, "Invalid table name '%s'\n", name);
        printout(prtbuf);
        return 0;
    }
    else
    {
        H5Gclose(project);
        return retcode;
    }

    return 1;
}

int printtable(char *name, char *chrname, hid_t project, unsigned int n0, unsigned int n1, int prj, int stride, FILE *outf)
{
    char buf[1001];
    hid_t dataset, datatype, dataspace;
    hsize_t dims[2], icount, ioff;
    int i, chr;
    int *valint;
    char **valstr;
    herr_t ret;

    icount = (hsize_t) (n1 - n0 + 1);
    ioff = (hsize_t) n0;

    //fprintf(outs, "%s %s %d %u %u %d %d\n", name, chrname, project, n0, n1, prj, stride);

    if (strcmp(name, "taxa") == 0)
    {
        dataset = H5Dopen(project, name, H5P_DEFAULT);
        if (dataset < 0)
        {
            sprintf(prtbuf, "Cannot open dataset '%s'\n", name);
            printout(prtbuf);
            return 0;
        }
        dataspace = H5Dget_space(dataset);
        if (dataspace < 0)
        {
            sprintf(buf, "Cannot open dataspace in dataset '%s'\n", name);
            printout(buf);
            H5Dclose(dataset);
            return 0;
        }
        //get dim
        ret = H5Sget_simple_extent_dims(dataspace, dims, NULL);
        if (dims[0] < n1)
        {
            sprintf(buf, "Ending index too large %u, maximum %llu\n", n1, dims[0]);
            printout(buf);
            H5Dclose(dataset);
            H5Sclose(dataspace);
            return 0;
        }
        datatype = H5Dget_type(dataset);
        if (datatype < 0)
        {
            sprintf(buf, "Cannot get datatype of dataset '%s'\n", name);
            H5Dclose(dataset);
            H5Sclose(dataspace);
            printout(buf);
            return 0;
        }
        if (read_str(ioff, icount, dataset, datatype, dataspace, &valstr) == -1)
        {
            H5Dclose(dataset);
            H5Sclose(dataspace);
            return -1;
        }
        for (i = 0; i < (unsigned int) icount; i += stride)
        {
            fprintf(outf, "%s\n", valstr[i]);
        }
        H5Dclose(dataset);
        H5Sclose(dataspace);
        freestrtable(valstr, icount);
    }
    else if (strcmp(name, "positions") == 0 || strcmp(name, "markers") == 0 || strcmp(name, "alleles") == 0)
    {
        chr = find_chr(chrname, prj);
        if (chr < 0)
        {
            printout("Invalid chromosome name\n");
            return 0;
        }
        sprintf(buf, "chr%d", chrinfo[chr].num);
        chr = H5Gopen(project, buf, H5P_DEFAULT);
        if (chr < 0)
        {
            printout("Cannot open chromosome\n");
            return 0;
        }
        if (strcmp(name, "markers") == 0)
        {
            if (!H5Lexists(chr, name, H5P_DEFAULT))
            {
                printout("Cannot open table markers\n");
                return 0;
            }
        }

        dataset = H5Dopen(chr, name, H5P_DEFAULT);
        if (dataset < 0)
        {
            sprintf(prtbuf, "Cannot open dataset '%s'\n", name);
            printout(prtbuf);
            H5Gclose(chr);
            return 0;
        }
        dataspace = H5Dget_space(dataset);
        if (dataspace < 0)
        {
            sprintf(buf, "Cannot open dataspace in dataset '%s'\n", name);
            printout(buf);
            H5Gclose(chr);
            H5Dclose(dataset);
            return 0;
        }
        ret = H5Sget_simple_extent_dims(dataspace, dims, NULL);
        if (dims[0] < n1)
        {
            sprintf(buf, "Ending index too large %u, maximum %llu\n", n1, dims[0]);
            printout(buf);
            H5Gclose(chr);
            H5Dclose(dataset);
            H5Sclose(dataspace);
            return 0;
        }
        datatype = H5Dget_type(dataset);
        if (datatype < 0)
        {
            sprintf(buf, "Cannot get datatype of dataset '%s'\n", name);
            H5Gclose(chr);
            H5Dclose(dataset);
            H5Sclose(dataspace);
            printout(buf);
            return 0;
        }
        if (strcmp(name, "positions") == 0)
        {
            if (read_int(ioff, icount, dataset, datatype, dataspace, &valint) == -1)
            {
                H5Gclose(chr);
                H5Dclose(dataset);
                H5Sclose(dataspace);
                return -1;
            }
            for (i = 0; i < (unsigned int) icount; i += stride)
            {
                fprintf(outf, "%u\n", valint[i]);
            }
            free(valint);
        }
        else //markers
        {
            if (read_str(ioff, icount, dataset, datatype, dataspace, &valstr) == -1)
            {
                H5Gclose(chr);
                H5Dclose(dataset);
                H5Sclose(dataspace);
                return -1;
            }
            for (i = 0; i < (unsigned int) icount; i += stride)
            {
                fprintf(outf, "%s\n", valstr[i]);
            }
            freestrtable(valstr, icount);
        }
        H5Gclose(chr);
        H5Dclose(dataset);
        H5Sclose(dataspace);
    }
    else
    {
        return -1;
    }
    return 1;
}

int mount()
{
    char buf[1001], name[256], pass[100];
    hid_t newfile, attr, project, type;
    int n = 0, i, j, ntsts, nn, order = -1, newmnt = 0;
    char buf1[257], orderstr[20];

    if (readline(pass, 100, __LINE__, 0) < 0)return -1;
    if (readline(name, 255, __LINE__, 1) < 0)return -1;
    if (readline(buf1, 256, __LINE__, 0) != -1)
    {
        printout("invalid input, expected end-of-command here\n");
        return -1;
    }
    if (check_pass(pass) == 0)
    {
        printout("invalid password\n");
        return 0;
    }
    //open file and analyze
    fprintf(flog, "Opening file %s for analyzing\n", name);
    fflush(flog);
    sprintf(prtbuf, "%s%s", hdf5dir, name);
    newfile = H5Fopen(prtbuf, H5P_DEFAULT, H5P_DEFAULT);
    if (newfile < 0)
    {
        sprintf(prtbuf, "Cannot open file %s\n", name);
        printout(prtbuf);
        return 0;
    }

    while (1)
    {
        sprintf(buf, "project%d", n);
        if (!H5Lexists(newfile, buf, H5P_DEFAULT))break;
        project = H5Gopen(newfile, buf, H5P_DEFAULT);
        fprintf(flog, "%s start\n", buf);
        fflush(flog);
        if (read_str_attr(project, "name", buf) < 0)
        {
            sprintf(buf, "Cannot read attribute 'name' of project%d in file %s\n", n, name);
            printout(buf);
            H5Gclose(project);
            H5Fclose(newfile);
            return 0;
        }
        fprintf(flog, "= %s\n", buf);
        fflush(flog);
        i = find_project(buf);
        fprintf(flog, "= %d\n", i);
        fflush(flog);
        if (i>-1)
        {
            sprintf(buf, "Project '%s' already exists in mounted file '%s'\n", projects[i].name, files[projects[i].mnt].name);
            printout(buf);
            H5Gclose(project);
            H5Fclose(newfile);
            return 0;
        }
        if (read_str_attr(project, "encoding", buf) < 0)
        {
            sprintf(buf, "Cannot read attribute 'encoding' of project%d in file %s\n", n, name);
            printout(buf);
            H5Gclose(project);
            H5Fclose(newfile);
            return 0;
        }
        if (read_str_attr(project, "genomeversion", prtbuf) < 0)
        {
            sprintf(buf, "Cannot read attribute 'genomeversion' of project%d in file %s\n", n, name);
            printout(buf);
            H5Gclose(project);
            H5Fclose(newfile);
            return 0;
        }
        if (read_str_attr(project, "orientation", orderstr) < 0)
        {
            sprintf(buf, "Cannot read attribute 'orientation' of project%d in file %s\n", n, name);
            printout(buf);
            H5Gclose(project);
            H5Fclose(newfile);
            return 0;
        }
        if (strcmp("pf", orderstr) == 0)
        {
            order = 1;
        }
        else if (strcmp("tf", orderstr) == 0)
        {
            order = 2;
        }
        else if (strcmp("both", orderstr) == 0)
        {
            order = 3;
        }
        else
        {
            sprintf(buf, "Invalid attribute 'orientation'='%s' of project%d in file %s\n", orderstr, n, name);
            printout(buf);
            H5Gclose(project);
            H5Fclose(newfile);
            return 0;
        }
        if (get_enc_size(buf) == -1)
        {
            sprintf(prtbuf, "Invalid attribute 'encoding'='%s' of project%d in file %s\n", buf, n, name);
            printout(prtbuf);
            H5Gclose(project);
            H5Fclose(newfile);
            return 0;
        }
        if (!H5Lexists(project, "taxa", H5P_DEFAULT))
        {
            sprintf(buf, "Missing 'taxa' table of project%d in file %s\n", n, name);
            printout(buf);
            H5Gclose(project);
            H5Fclose(newfile);
            return 0;
        }
        if (!H5Lexists(project, "taxa_sv", H5P_DEFAULT) || !H5Lexists(project, "taxa_si", H5P_DEFAULT))
        {
            sprintf(buf, "WARNING: 'taxa' table of project%d in file %s is not indexed - text searches are disabled\n", n, name);
            printout(buf);
        }
        if (read_chrinfo(project, -1, -1, 1, 1, order) < 0)
        {
            sprintf(prtbuf, "Error reading project no %d chromosomes \n", n);
            printout(prtbuf);
            H5Gclose(project);
            H5Fclose(newfile);
            return 0;
        }
        H5Gclose(project);
        n++;
        fprintf(flog, "%s end OK\n", buf);
        fflush(flog);
    }
    if (n == 0)
    {
        sprintf(buf, "No projects found in the file %s\n", name);
        printout(buf);
        H5Fclose(newfile);
        return 0;
    }
    fprintf(flog, "Verified %d projects, all OK\nAdding file\n", n);
    fflush(flog);

    //add file
    if (nfiles > 0)
    {
        newmnt = -1;
        for (i = 0; i < nfiles; i++)
        {
            ntsts = -1;
            for (j = 0; j < nfiles; j++)
            {
                if (files[j].mnt == i)ntsts = 1;
            }
            if (ntsts == -1)
            {
                newmnt = i;
                break;
            }
        }
        if (newmnt == -1)newmnt = nfiles;
        nfiles++;
        files = realloc(files, nfiles * sizeof (files_type));
        filesh = realloc(filesh, nfiles * sizeof (hid_t));
    }
    else
    {
        nfiles = 1;
        newmnt = 0;
    }
    strcpy(files[nfiles - 1].name, name);
    files[nfiles - 1].mnt = newmnt;
    filesh[nfiles - 1] = newfile;

    fprintf(flog, "Adding projects\n");
    fflush(flog);
    //add projects
    nn = nprojects;
    nprojects += n;
    if (nn == 0)
    {
        nprojects = n;
    }
    else
    {
        projects = realloc(projects, nprojects * sizeof (projects_type));
    }
    for (i = 0; i < n; i++)
    {
        sprintf(buf, "/project%d", i);
        project = H5Gopen(newfile, buf, H5P_DEFAULT);
        attr = H5Aopen_by_name(project, ".", "name", H5P_DEFAULT, H5P_DEFAULT);
        type = H5Tcopy(H5T_C_S1);
        H5Tset_size(type, 256);
        H5Aread(attr, type, buf);
        projects[nn + i].mnt = newmnt;
        projects[nn + i].num = i;
        strcpy(projects[nn + i].name, buf);
        read_str_attr(project, "encoding", buf);
        projects[nn + i].enc = get_enc_size(buf);
        strcpy(projects[nn + i].encoding, buf);
        read_str_attr(project, "orientation", buf);
        if (strcmp("pf", buf) == 0)
        {
            order = 1;
        }
        else if (strcmp("tf", buf) == 0)
        {
            order = 2;
        }
        else
        {
            order = 3;
        }
        projects[nn + i].order = read_chrinfo(project, newmnt, nn + i, 0, 0, order);
        H5Gclose(project);
        H5Aclose(attr);
    }
    //write files and projects to file
    write_root_props('F', 0);
    write_root_props('P', 0);

    sprintf(buf, "/mnt%d", nfiles - 1);
    if (!H5Lexists(file, buf, H5P_DEFAULT))
    {
        H5Gclose(H5Gcreate(file, buf, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
    }
    H5Fflush(file, H5F_SCOPE_LOCAL);
    if (H5Fmount(file, buf, newfile, H5P_DEFAULT) < 0)
    {
        sprintf(prtbuf, "Cannot mount file %s on group %s", name, buf);
        error(prtbuf, __LINE__);
    }
    fflush(outs);
    return 1;
}

int delproject(int n)
{
    int i, nn, j;
    projects_type *projects1;
    int nprojects1;
    int nchrinfo1;
    chrinfo_type *chrinfo1;

    //clean up projects
    nn = 1;
    if (nprojects - nn > 1)
    {
        nprojects1 = nprojects - nn;
        projects1 = (projects_type *) malloc(nprojects1 * sizeof (projects_type));
        j = 0;
        for (i = 0; i < nprojects; i++)
        {
            if (i != n)
            {
                memcpy(&projects1[j], &projects[i], sizeof (projects_type));
                j++;
            }
        }
        free(projects);
        projects = projects1;
        nprojects = nprojects1;
    }
    else
    {
        if (nprojects > 1)projects = realloc(projects, sizeof (projects_type));
        nprojects = 0;
        projects[0].mnt = -1;
        projects[0].num = -1;
        projects[0].name[0] = '\0';
    }
    //clean up chrinfo
    nn = 0;
    for (i = 0; i < nchrinfo; i++)
        if (chrinfo[i].prj == n)nn++;
    if (nprojects - nn > 1)
    {
        nchrinfo1 = nchrinfo - nn;
        chrinfo1 = (chrinfo_type *) malloc(nchrinfo1 * sizeof (chrinfo_type));
        j = 0;
        for (i = 0; i < nchrinfo; i++)
        {
            if (chrinfo[i].prj != n)
            {
                memcpy(&chrinfo1[j], &chrinfo[i], sizeof (chrinfo_type));
                j++;
            }
        }
        free(chrinfo);
        chrinfo = chrinfo1;
        nchrinfo = nchrinfo1;
    }
    else
    {
        if (nchrinfo > 0)free(chrinfo);
        nchrinfo = 0;
    }
    //write files and projects to file
    write_root_props('P', 1);
    H5Fflush(file, H5F_SCOPE_LOCAL);
    return 1;
}

int umount1(int n)
{
    int i, nn, j;
    projects_type *projects1;
    int nprojects1;
    files_type *files1;
    int nfiles1;
    hid_t *filesh1;
    int nchrinfo1;
    chrinfo_type *chrinfo1;

    //n = find_file(name);
    //unmount the file
    //clean up files filesh
    if (nfiles > 1)
    {
        nfiles1 = nfiles - 1;
        files1 = (files_type *) malloc(nfiles1 * sizeof (files_type));
        filesh1 = (hid_t *) malloc(nfiles1 * sizeof (hid_t));
        j = 0;
        for (i = 0; i < nfiles; i++)
        {
            if (i != n)
            {
                memcpy(&files1[j], &files[i], sizeof (files_type));
                filesh1[j] = filesh[i];
                j++;
            }
        }
        free(files);
        files = files1;
        free(filesh);
        filesh = filesh1;
        nfiles = nfiles1;
    }
    else
    {
        nfiles = 0;
        files[0].mnt = -1;
        files[0].name[0] = '\0';
    }
    //clean up projects
    nn = 0;
    for (i = 0; i < nprojects; i++)
        if (projects[i].mnt == n)nn++;
    if (nprojects - nn > 1)
    {
        nprojects1 = nprojects - nn;
        projects1 = (projects_type *) malloc(nprojects1 * sizeof (projects_type));
        j = 0;
        for (i = 0; i < nprojects; i++)
        {
            if (projects[i].mnt != n)
            {
                memcpy(&projects1[j], &projects[i], sizeof (projects_type));
                j++;
            }
        }
        free(projects);
        projects = projects1;
        nprojects = nprojects1;
    }
    else
    {
        if (nprojects > 1)projects = realloc(projects, sizeof (projects_type));
        nprojects = 0;
        projects[0].mnt = -1;
        projects[0].num = -1;
        projects[0].name[0] = '\0';
    }
    //clean up chrinfo
    nn = 0;
    for (i = 0; i < nchrinfo; i++)
        if (chrinfo[i].mnt == n)nn++;
    if (nprojects - nn > 1)
    {
        nchrinfo1 = nchrinfo - nn;
        chrinfo1 = (chrinfo_type *) malloc(nchrinfo1 * sizeof (chrinfo_type));
        j = 0;
        for (i = 0; i < nchrinfo; i++)
        {
            if (chrinfo[i].mnt != n)
            {
                memcpy(&chrinfo1[j], &chrinfo[i], sizeof (chrinfo_type));
                j++;
            }
        }
        free(chrinfo);
        chrinfo = chrinfo1;
        nchrinfo = nchrinfo1;
    }
    else
    {
        if (nchrinfo > 0)free(chrinfo);
        nchrinfo = 0;
    }
    //write files and projects to file
    write_root_props('F', 1);
    write_root_props('P', 1);
    H5Fflush(file, H5F_SCOPE_LOCAL);
    return 1;
}

int umount()
{
    char buf[1001], name[256], pass[100];
    int n = 0, i, nn, j;
    projects_type *projects1;
    int nprojects1;
    files_type *files1;
    int nfiles1;
    hid_t *filesh1;
    int nchrinfo1;
    chrinfo_type *chrinfo1;

    if (readline(pass, 100, __LINE__, 0) < 0)return -1;
    if (readline(name, 255, __LINE__, 1) < 0)return -1;
    if (readline(buf, 256, __LINE__, 0) != -1)
    {
        printout("invalid input, expected end-of-command here\n");
        return -1;
    }
    if (check_pass(pass) == 0)
    {
        printout("invalid password\n");
        return 0;
    }
    n = find_file(name);
    if (n < 0)
    {
        printout("Invalid file name\n");
        return 0;
    }
    //unmount the file
    sprintf(buf, "/mnt%d", files[n].mnt);
    H5Funmount(file, buf);
    H5Fclose(filesh[n]);
    //clean up files filesh
    if (nfiles > 1)
    {
        nfiles1 = nfiles - 1;
        files1 = (files_type *) malloc(nfiles1 * sizeof (files_type));
        filesh1 = (hid_t *) malloc(nfiles1 * sizeof (hid_t));
        j = 0;
        for (i = 0; i < nfiles; i++)
        {
            if (i != n)
            {
                memcpy(&files1[j], &files[i], sizeof (files_type));
                files1[j].mnt = j;
                filesh1[j] = filesh[i];
                if (i > n)
                {
                    sprintf(buf, "/mnt%d", i);
                    H5Funmount(file, buf);
                    sprintf(buf, "/mnt%d", j);
                    if (H5Fmount(file, buf, filesh1[j], H5P_DEFAULT) < 0)
                    {
                        sprintf(prtbuf, "Cannot mount file %s on group %s", name, buf);
                        error(prtbuf, __LINE__);
                    }
                }
                j++;
            }
        }
        free(files);
        files = files1;
        free(filesh);
        filesh = filesh1;
        nfiles = nfiles1;
    }
    else
    {
        nfiles = 0;
        files[0].mnt = -1;
        files[0].name[0] = '\0';
    }
    //clean up projects
    nn = 0;
    for (i = 0; i < nprojects; i++)
        if (projects[i].mnt == n)nn++;
    if (nprojects - nn >= 1)
    {
        nprojects1 = nprojects - nn;
        projects1 = (projects_type *) malloc(nprojects1 * sizeof (projects_type));
        j = 0;
        for (i = 0; i < nprojects; i++)
        {
            if (projects[i].mnt != n)
            {
                memcpy(&projects1[j], &projects[i], sizeof (projects_type));
                if (projects[i].mnt > n)projects1[j].mnt = projects[i].mnt - 1;
                j++;
            }
        }
        free(projects);
        projects = projects1;
        nprojects = nprojects1;
    }
    else
    {
        if (nprojects > 1)projects = realloc(projects, sizeof (projects_type));
        nprojects = 0;
        projects[0].mnt = -1;
        projects[0].num = -1;
        projects[0].name[0] = '\0';
    }
    //clean up chrinfo
    nn = 0;
    for (i = 0; i < nchrinfo; i++)
        if (chrinfo[i].mnt == n)nn++;
    if (nchrinfo - nn >= 1)
    {
        nchrinfo1 = nchrinfo - nn;
        chrinfo1 = (chrinfo_type *) malloc(nchrinfo1 * sizeof (chrinfo_type));
        j = 0;
        for (i = 0; i < nchrinfo; i++)
        {
            if (chrinfo[i].mnt != n)
            {
                memcpy(&chrinfo1[j], &chrinfo[i], sizeof (chrinfo_type));
                if (chrinfo[i].mnt > n)chrinfo1[j].mnt = chrinfo[i].mnt - 1;
                j++;
            }
        }
        free(chrinfo);
        chrinfo = chrinfo1;
        nchrinfo = nchrinfo1;
    }
    else
    {
        if (nchrinfo > 0)free(chrinfo);
        nchrinfo = 0;
    }
    //write files and projects to file
    write_root_props('F', 1);
    write_root_props('P', 1);
    H5Fflush(file, H5F_SCOPE_LOCAL);
    fflush(outs);
    return 1;
}

int useradd()
{
    char buf[1001], name[256], pass[100], upass[100];

    if (readline(pass, 100, __LINE__, 0) < 0)return -1;
    if (readline(name, 255, __LINE__, 1) < 0)return -1;
    if (readline(upass, 255, __LINE__, 0) < 0)return -1;
    if (readline(buf, 256, __LINE__, 0) != -1)
    {
        printout("invalid input, expected end-of-command here\n");
        return -1;
    }
    if (check_pass(pass) == 0)
    {
        printout("invalid password\n");
        return 0;
    }

    if (nusers > 0)
    {
        nusers++;
        users = realloc(users, nusers * sizeof (userdata_type));
        accesstable = realloc(accesstable, nusers * sizeof (char *));
    }
    else
    {
        nusers = 1;
    }

    strcpy(users[nusers - 1].username, name);
    set_user_pass(nusers - 1, upass);
    accesstable[nusers - 1] = NULL;
    write_root_props('U', 0);
    write_root_props('A', 0);
    fprintf(outs, "%d\n", nusers);

    return 1;
}

int userdel()
{
    char buf[1001], name[256], pass[100];
    int n = 0, j, i, nusers1;
    userdata_type *users1;
    char **access1;


    if (readline(pass, 100, __LINE__, 0) < 0)return -1;
    if (readline(name, 255, __LINE__, 1) < 0)return -1;
    if (readline(buf, 256, __LINE__, 0) != -1)
    {
        printout("invalid input, expected end-of-command here\n");
        return -1;
    }
    if (check_pass(pass) == 0)
    {
        printout("invalid password\n");
        return 0;
    }

    n = find_user(name);
    if (n == -1)
    {
        printout("invalid user name\n");
        return 0;
    }

    if (nusers > 1)
    {
        nusers1 = nusers - 1;
        users1 = (userdata_type *) malloc(nusers1 * sizeof (userdata_type));
        access1 = (char **) malloc(nusers1 * sizeof (char *));
        j = 0;
        for (i = 0; i < nusers; i++)
        {
            if (i != n)
            {
                memcpy(&users1[j], &users[i], sizeof (userdata_type));
                access1[j] = accesstable[i];
                j++;
            }
        }
        free(users);
        freestrtable(accesstable, nusers);
        users = users1;
        accesstable = access1;
        nusers = nusers1;
    }
    else
    {
        nusers = 0;
        memset(users[0].username, '\0', 50);
        accesstable[0] = NULL;
    }

    write_root_props('U', 1);
    write_root_props('A', 1);

    return 1;
}

int userpass()
{
    char buf[1001], name[256], pass[100], newpass[100];

    if (readline(name, 255, __LINE__, 1) < 0)return -1;
    if (readline(pass, 100, __LINE__, 0) < 0)return -1;
    if (readline(newpass, 100, __LINE__, 0) < 0)return -1;
    if (readline(buf, 256, __LINE__, 0) != -1)
    {
        printout("invalid input, expected end-of-command here\n");
        return -1;
    }
    int usr = check_user_pass(name, pass);
    if (usr == -1)
    {
        //check if it is admin pass
        if (check_pass(pass) == 0)
        {
            printout("invalid current password\n");
            return 0;
        }
        usr = get_user(name);
    }

    set_user_pass(usr, newpass);
    write_root_props('U', 0);

    return 1;
}

int useracc()
{
    char buf[1001], name[256], pass[100];
    int n = 0, len = 0, i = 0, usr;
    char *acc = NULL;

    if (readline(pass, 100, __LINE__, 0) < 0)return -1;
    if (readline(name, 255, __LINE__, 1) < 0)return -1;
    while (readline(buf, 256, __LINE__, 0) != -1)
    {
        chomp(buf);
        if (n == 0)
        {
            acc = malloc(strlen(buf) + 3);
            sprintf(acc, "\n%s\n", buf);
            len = strlen(buf) + 3;
            i = len - 1;
        }
        else
        {
            len = len + strlen(buf) + 2;
            acc = realloc(acc, len);
            sprintf(&acc[i], "\n%s\n", buf);
            i = len - 1;
        }
        n++;
    }
    if (check_pass(pass) == 0)
    {
        printout("invalid password\n");
        return 0;
    }
    usr = find_user(name);
    if (usr == -1)
    {
        printout("invalid user name\n");
        return 0;
    }
    if (accesstable[usr] != NULL)free(accesstable[usr]);
    accesstable[usr] = acc;
    write_root_props('A', 1);
    return 1;
}

int pinfo()
{
    int i, n, order = 0;
    hsize_t dims[2];
    char buf[1000];
    char buf1[257];
    hid_t project, chromosome, dataset, filespace;
    herr_t ret;

    if (readline(buf, 256, __LINE__, 1) < 0)return -1;
    if (readline(buf1, 256, __LINE__, 0) != -1)
    {
        printout("invalid input, expected end-of-command here\n");
        return -1;
    }
    n = find_project(buf);
    if (n == -1)
    {
        printout("Invalid project name\n");
        return 0;
    }
    sprintf(buf, "/mnt%d/project%d\t%s\n", projects[n].mnt, projects[n].num, projects[n].name);
    printout(buf);

    sprintf(buf, "/mnt%d/project%d", projects[n].mnt, projects[n].num);
    project = H5Gopen(file, buf, H5P_DEFAULT);
    if (project < 0)error("Cannnot open project", __LINE__);

    if (read_str_attr(project, "type", buf) < 0)
    {
        fprintf(outs, "cannot read attribute 'type'\n");
    }
    else
    {
        fprintf(outs, "type\t%s\n", buf);
    }
    if (read_str_attr(project, "version", buf) < 0)
    {
        fprintf(outs, "cannot read attribute 'version'\n");
    }
    else
    {
        fprintf(outs, "version\t%s\n", buf);
    }
    if (read_str_attr(project, "genomeversion", buf) < 0)
    {
        fprintf(outs, "cannot read attribute 'genomeversion'\n");
    }
    else
    {
        fprintf(outs, "genomeversion\t%s\n", buf);
    }
    if (read_str_attr(project, "orientation", buf) < 0)
    {
        fprintf(outs, "cannot read attribute 'orientation'\n");
    }
    else
    {
        fprintf(outs, "orientation\t%s\n", buf);
        order = 1;
        if (strcmp(buf, "both") == 0)order = 3;
        if (strcmp(buf, "tf") == 0)order = 2;
    }
    if (read_int_attr(project, "status", &i) < 0)
    {
        fprintf(outs, "cannot read attribute 'status'\n");
    }
    else
    {
        fprintf(outs, "status\t%d\n", i);
    }
    if (read_str_attr(project, "coordinate_system", buf) < 0)
    {
        fprintf(outs, "cannot read attribute 'coordinate_system'\n");
    }
    else
    {
        fprintf(outs, "coordinate_system\t%s\n", buf);
    }
    if (read_str_attr(project, "encoding", buf) < 0)
    {
        fprintf(outs, "cannot read attribute 'encoding'\n");
    }
    else
    {
        fprintf(outs, "encoding\t%s\n", buf);
    }
    dataset = H5Dopen(project, "taxa", H5P_DEFAULT);
    if (dataset < 0)
    {
        fprintf(outs, "=\tCannot open dataset taxa\n");
    }
    else
    {
        filespace = H5Dget_space(dataset);
        ret = H5Sget_simple_extent_dims(filespace, dims, NULL);
        fprintf(outs, "=\ttaxa\t%lli\t", dims[0]);
        H5Sclose(filespace);
        H5Dclose(dataset);
        if (!H5Lexists(project, "taxa_sv", H5P_DEFAULT) || !H5Lexists(project, "taxa_si", H5P_DEFAULT))
        {
            fprintf(outs, "NOT indexed\n");
        }
        else
        {
            fprintf(outs, "indexed\n");
        }
    }
    //read chromosomes
    i = 0;
    while (1)
    {
        i++;
        sprintf(buf, "chr%d", i);
        if (!H5Lexists(project, buf, H5P_DEFAULT))break;
        chromosome = H5Gopen(project, buf, H5P_DEFAULT);
        if (read_str_attr(chromosome, "name", buf) < 0)
        {
            fprintf(outs, "chromosome\t%d\tchr%d\n", i, i);
        }
        else
        {
            fprintf(outs, "chromosome\t%d\t%s\n", i, buf);
        }
        if (order == 1 || order == 3)
        {
            sprintf(buf, "%s_pf", DATAARRAY);
            dataset = H5Dopen(chromosome, buf, H5P_DEFAULT);
            if (dataset < 0)
            {
                sprintf(prtbuf, "=\tCannot open dataset '%s'\n", buf);
                fprintf(outs, prtbuf);
            }
            else
            {
                filespace = H5Dget_space(dataset);
                ret = H5Sget_simple_extent_dims(filespace, dims, NULL);
                fprintf(outs, "=\t%s\t%lli\t%lli\n", buf, dims[0], dims[1]);
                H5Sclose(filespace);
                H5Dclose(dataset);
            }
        }
        if (order == 2 || order == 3)
        {
            sprintf(buf, "%s_tf", DATAARRAY);
            dataset = H5Dopen(chromosome, buf, H5P_DEFAULT);
            if (dataset < 0)
            {
                sprintf(prtbuf, "=\tCannot open dataset '%s'\n", buf);
                fprintf(outs, prtbuf);
            }
            else
            {
                filespace = H5Dget_space(dataset);
                ret = H5Sget_simple_extent_dims(filespace, dims, NULL);
                fprintf(outs, "=\t%s\t%lli\t%lli\n", buf, dims[0], dims[1]);
                H5Sclose(filespace);
                H5Dclose(dataset);
            }
        }
        dataset = H5Dopen(chromosome, "positions", H5P_DEFAULT);
        if (dataset < 0)
        {
            fprintf(outs, "=\tCannot open dataset 'positions'\n");
        }
        else
        {
            filespace = H5Dget_space(dataset);
            ret = H5Sget_simple_extent_dims(filespace, dims, NULL);
            fprintf(outs, "=\tpositions\t%lli\n", dims[0]);
            H5Sclose(filespace);
            H5Dclose(dataset);
        }

        dataset = H5Dopen(chromosome, "alleles", H5P_DEFAULT);
        if (dataset < 0)
        {
            fprintf(outs, "=\tCannot open dataset 'alleles'\n");
        }
        else
        {
            filespace = H5Dget_space(dataset);
            ret = H5Sget_simple_extent_dims(filespace, dims, NULL);
            fprintf(outs, "=\talleles\t%lli\n", dims[0]);
            H5Sclose(filespace);
            H5Dclose(dataset);
        }

        if (!H5Lexists(chromosome, "markers", H5P_DEFAULT))
        {
            fprintf(outs, "=\tDataset 'markers' not present\n");
        }
        else
        {
            dataset = H5Dopen(chromosome, "markers", H5P_DEFAULT);
            filespace = H5Dget_space(dataset);
            ret = H5Sget_simple_extent_dims(filespace, dims, NULL);
            fprintf(outs, "=\tmarkers\t%lli\t", dims[0]);
            H5Sclose(filespace);
            H5Dclose(dataset);
            if (!H5Lexists(chromosome, "markers_sv", H5P_DEFAULT) || !H5Lexists(chromosome, "markers_si", H5P_DEFAULT))
            {
                fprintf(outs, "NOT indexed\n");
            }
            else
            {
                fprintf(outs, "indexed\n");
            }
        }
        H5Gclose(chromosome);
    }
    fflush(outs);
    return 1;
}

int plist()
{
    int i;
    char buf[257];

    if (readline(buf, 256, __LINE__, 0) != -1)
    {
        printout("invalid input, expected end-of-command here\n");
        return -1;
    }

    fprintf(outs, "%d projects attached\n", nprojects);
    for (i = 0; i < nprojects; i++)
    {
        fprintf(outs, "/mnt%d/project%d\t%s\n", projects[i].mnt, projects[i].num, projects[i].name);
    }
    fflush(outs);
    return 1;
}

int flist()
{
    int i;
    char buf[257];

    if (readline(buf, 256, __LINE__, 0) != -1)
    {
        printout("invalid input, expected end-of-command here\n");
        return -1;
    }

    fprintf(outs, "%d files attached\n", nfiles);
    for (i = 0; i < nfiles; i++)
    {
        fprintf(outs, "/mnt%d\t%s\n", i, files[i].name);
    }
    fflush(outs);
    return 1;
}

int login()
{
    char buf[1001], name[256], pass[100];

    if (readline(name, 255, __LINE__, 1) < 0)return -1;
    if (readline(pass, 100, __LINE__, 0) < 0)return -1;
    if (readline(buf, 256, __LINE__, 0) != -1)
    {
        printout("invalid input, expected end-of-command here\n");
        return -1;
    }
    int usr = check_user_pass(name, pass);
    if (usr == -1)
    {
        if (check_pass(pass) != 0 && strcmp(name, "serveradmin") == 0)
        {
            printout("OK\n");
        }
        else
        {
            printout("ERR\n");
        }
    }
    else
    {
        printout("OK\n");
        printout(&accesstable[usr][1]);
    }
    return 1;
}

int aflist()
{
    int i;
    char buf[257], pass[100];
    FILE *in;

    if (readline(pass, 100, __LINE__, 0) < 0)return -1;
    if (readline(buf, 256, __LINE__, 0) != -1)
    {
        printout("invalid input, expected end-of-command here\n");
        return -1;
    }

    if (check_pass(pass) == 0)
    {
        printout("invalid password\n");
        return 0;
    }

    sprintf(prtbuf, "ls %s/data > %s/ooo", hdf5conf, tmpdir);
    i = system(prtbuf);
    sprintf(prtbuf, "%s/ooo", tmpdir);
    if ((in = fopen(prtbuf, "r")) == NULL)
    {
        printout("cannot list files\n");
        return 0;
    }
    while (fgets(prtbuf, 4999, in))
    {
        chomp(prtbuf);
        if (strcmp(prtbuf, "hdf5serverroot.h5") != 0 && strcmp(prtbuf, "hdf5pipe_out") != 0 && strcmp(prtbuf, "hdf5pipe_in") != 0 && strcmp(prtbuf, "hdf5.lock") != 0)
        {
            sprintf(buf, "%s\n", prtbuf);
            printout(buf);
        }
        //sprintf(buf, "*%s*\n", prtbuf);
        //printout(buf);
    }

    fclose(in);
    sprintf(prtbuf, "rm -f %s/ooo", tmpdir);
    i = system(prtbuf);
    fflush(outs);
    return 1;
}

int ulist()
{
    int i;
    char buf[257], pass[100];

    if (readline(pass, 100, __LINE__, 0) < 0)return -1;
    if (readline(buf, 256, __LINE__, 0) != -1)
    {
        printout("invalid input, expected end-of-command here\n");
        return -1;
    }

    if (check_pass(pass) == 0)
    {
        printout("invalid password\n");
        return 0;
    }

    if (nusers == 0)return 1;
    for (i = 0; i < nusers; i++)
    {
        sprintf(prtbuf, "USER %s %s\n", users[i].username, users[i].password);
        printout(prtbuf);
        if (accesstable[i] == NULL)
        {
            //do nothing, no access to projects
        }
        else
        {
            printout(&accesstable[i][1]);
        }
    }

    return 1;
}

int finfo()
{
    char buf[256];
    char buf1[257];
    int i, n;

    if (readline(buf, 256, __LINE__, 1) < 0)return -1;
    if (readline(buf1, 256, __LINE__, 0) != -1)
    {
        printout("invalid input, expected end-of-command here\n");
        return -1;
    }
    n = find_file(buf);
    if (n == -1)
    {
        printout("Invalid file name\n");
        return 0;
    }
    sprintf(buf, "/mnt%d\t%s\n", n, files[n].name);
    printout(buf);
    for (i = 0; i < nprojects; i++)
    {
        if (projects[i].mnt == files[n].mnt)
        {
            fprintf(outs, "/mnt%d/project%d\t%s\n", projects[i].mnt, projects[i].num, projects[i].name);
        }
    }
    fflush(outs);
    return 1;
}

void shutdown()
{
    int i;
    char buf[100];

    for (i = 0; i < nfiles; i++)
    {
        sprintf(buf, "/mnt%d", files[i].mnt);
        H5Funmount(file, buf);
        H5Fclose(filesh[i]);
    }
    H5Fclose(file);
    if (comm_pid != -1)kill(comm_pid, 9);
}
