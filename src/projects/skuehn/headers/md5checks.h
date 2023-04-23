/**
   \file md5checks.cpp
   File which provides an easy function to check whether two files have equal MD5-sum or not. It allows for quick checks when/where files get corrupted. Whenever a pointer is returned the user has to free the memory again, otherwise there will be a leak. For this function to work parts of the open-ssl library have to be linked the following way: -lssl -lcrypto
   
   \todo Get it more clear and improve the memory management. Up to now the user might get a pointer and has to take care of freeing the memory himself.
   
  \author Stefan Kühn
  \date 01/07/2013

*/

#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <openssl/md5.h>



//Print the MD5 sum as hex-digits.
void print_md5_sum(unsigned char* md) 
{
    int i;
    for(i=0; i <MD5_DIGEST_LENGTH; i++) 
    {
            printf("%02x",md[i]);
    }
    printf("\n");
}

//Get the size of the file by its file descriptor
unsigned long get_size_by_fd(int fd) 
{
    struct stat statbuf;
    if(fstat(fd, &statbuf) < 0) exit(-1);
    return statbuf.st_size;
}

//Get the md5sum and return a pointer to it
unsigned char* get_md5sum (const char* filename)
{
    int file_descript;
    unsigned long file_size;
    char* file_buffer;
    unsigned char* result;
    result = (unsigned char*)malloc(MD5_DIGEST_LENGTH*sizeof(unsigned char));

    printf("using file:\t%s\n", filename);

    file_descript = open(filename, O_RDONLY);
    if(file_descript < 0) exit(-1);

    file_size = get_size_by_fd(file_descript);
    printf("file size:\t%lu\n", file_size);

    file_buffer = (char*)  mmap(0,file_size, PROT_READ,MAP_SHARED,file_descript,0);
    MD5((unsigned char*) file_buffer, file_size, result);
    munmap(file_buffer, file_size);
    
    return result;

  
}

//Compare two MD5-sums
bool md5_compare(unsigned char* md5sum1, unsigned char* md5sum2)
{
  if(strcmp((char*) md5sum1,(char*) md5sum2)==0)
    return true;
  else 
    return false;
}

//Routine which takes two filenames and checks of the MD5-sums of both files are equal.
bool md5_file_compare(char filename1[], char filename2[])
{
    unsigned char *md5sum1,*md5sum2;
    bool equal;
    
    md5sum1 = get_md5sum(filename1);
    md5sum2 = get_md5sum(filename2);
    
    equal = md5_compare(md5sum1,md5sum2);
    

    if(equal)
    {
      printf("MD5-sum file1:\t");
      print_md5_sum(md5sum1);
      printf("MD5-sum file2:\t");
      print_md5_sum(md5sum2);
      printf("Sums are equal\n");
    }
    else
    {
      printf("MD5-sum file1:\t");
      print_md5_sum(md5sum1);
      printf("MD5-sum file2:\t");
      print_md5_sum(md5sum2);
      printf("Sums are NOT equal\n");
    }      
    
    return md5_compare(md5sum1,md5sum2);
    
    free(md5sum1);
    free(md5sum2);
}

