#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

#include <my_timer.h>

#include <my_omp.h>

struct LinkedList_t
{
   LinkedList_t *next;

   int value;

   LinkedList_t() : next(NULL), value(0) {}
};

LinkedList_t * buildList (const int n)
{
   LinkedList_t *head = new LinkedList_t;
   head->value = 0;
   head->next = NULL;

   LinkedList_t *p = head;
   for (int i = 1; i < n; ++i, p = p->next)
   {
      LinkedList_t *item = new LinkedList_t;
      item->value = i;
      item->next  = NULL;
      p->next = item;
   }

   return head;
}

void destroyList (LinkedList_t *head)
{
   LinkedList_t *p = head;
   while (p)
   {
      LinkedList_t *next = p->next;
      delete p;
      p = next;
   }
}

void printList (LinkedList_t *head)
{
   for (LinkedList_t *p = head; p != NULL; p = p->next)
   {
      printf("item=%x value=%d\n", p, p->value);
   }
}

int sumList;

void doSomething (LinkedList_t *item)
{
   int old_value = item->value;

   item->value = 0;
   for (int i = 0; i < 10*old_value; ++i)
      item->value += sqrt(i);

   int threadId = 0;
#ifdef _OPENMP
   threadId = omp_get_thread_num();
#endif

   printf("item=%x value=%d %d %d\n", item, old_value, item->value, threadId);

   item->value = old_value;

   // Update the global counter one thread at a time.
   #pragma omp critical
   {
      sumList += item->value;
   }
}

void computeList (LinkedList_t *head)
{
   LinkedList_t *p = head;

   sumList = 0;

   for (LinkedList_t *p = head; p != NULL; p = p->next)
   {
      doSomething(p);
   }

   printf("computeList: sumList=%d\n", sumList);
}
void computeListOMP (LinkedList_t *head)
{
   LinkedList_t *p = head;

   sumList = 0;

   #pragma omp parallel default(shared) shared(p)
   {
      int threadId = 0;
#ifdef _OPENMP
      threadId = omp_get_thread_num();
#endif
      while (true)
      {
         LinkedList_t *myItem = NULL;
         #pragma omp critical
         {
            if (p)
            {
               myItem = p;
               p = p->next;
            }
         }

         if (myItem)
         {
            doSomething(myItem);
         }
         else
            break;
      }
   }

   printf("computeListOMP: sumList=%d\n", sumList);
}

#if defined(_OPENMP) && (_MY_OPENMP_MAJOR >= 3)
void computeListTask (LinkedList_t *head)
{
   sumList = 0;

   LinkedList_t *p = head;

   #pragma omp parallel default(shared)
   {
      #pragma omp single nowait
      {
         while(p)
         {
            #pragma omp task firstprivate(p)
            { doSomething(p); }
            p = p->next;
         }
      }
   }

   printf("computeListTask: sumList=%d\n", sumList);
}
#endif

int main(int argc, char *argv[])
{
   int num_threads = 1;

#ifdef _OPENMP
   fprintf(stderr,"OpenMP specification c. %s %d.%d\n", TOSTRING(_OPENMP), _MY_OPENMP_MAJOR, _MY_OPENMP_MINOR);
# ifdef _MY_OPENMP_HAS_ATOMICS
   fprintf(stderr,"OpenMP has atomics!\n");
# endif

   #pragma omp parallel default(shared)
   {
      #pragma omp master
      num_threads = omp_get_num_threads();
   }
   fprintf(stderr,"num_threads = %d\n", num_threads);
#endif

   int n = 10;

   for (int i = 1; i < argc; ++i)
   {
#define check_index(i,str) \
   if ((i) >= argc) \
      { fprintf(stderr,"Missing 2nd argument for %s\n", str); return 1; }

      if ( strcmp(argv[i],"-h") == 0 || strcmp(argv[i],"--help") == 0)
      {
         fprintf(stderr,"omp_list --help|-h --length|-n\n");
         return 1;
      }
      else if (strcmp(argv[i],"--length") == 0 || strcmp(argv[i],"-n") == 0)
      {
         check_index(i+1,"--length|-n");
         i++;
         if (isdigit(*argv[i]))
            n = atoi( argv[i] );
      }
   }

   printf("length: %d\n", n);

   LinkedList_t *head = buildList(n);

   //printList(head);

   computeList(head);

   computeListOMP(head);

#if defined(_OPENMP) && (_MY_OPENMP_MAJOR >= 3)
   computeListTask(head);
#endif

   destroyList(head);

   return 0;
}
