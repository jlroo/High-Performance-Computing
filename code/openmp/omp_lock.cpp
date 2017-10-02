#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

#include <my_omp.h>
#include <my_timer.h>

template <typename T>
struct List
{
   struct Node
   {
      Node *next, *prev;
      T item;
      int rank;

      Node() : next(NULL), prev(NULL), item(), rank(-1) {}
   };

   Node *head, *tail;

#ifdef _OPENMP
   omp_lock_t lock;
#endif

   List() : head(NULL), tail(NULL)
   {
      head = new Node;
      tail = new Node;
      head->next = tail;
      tail->prev = head;

#ifdef _OPENMP
      omp_init_lock( &this->lock );
#endif
   }

   ~List()
   {
      Node *p = head->next;
      while (p != tail)
      {
         Node *next = p->next;
         delete p;
         p = next;
      }

      delete head;
      delete tail;

#ifdef _OPENMP
      omp_destroy_lock( &this->lock );
#endif
   }

   bool isEmpty (void)
   {
      bool res;

#ifdef _OPENMP
      omp_set_lock( &lock );
#endif

      res = (head->next == tail);

#ifdef _OPENMP
      omp_unset_lock( &lock );
#endif

      return res;
   }

   T& front (void)
   {
      if (head->next != tail)
         return head->next->item;
   }

   T& back (void)
   {
      if (tail->prev != head)
         return tail->prev->item;
   }

   bool popFront (T *item)
   {
      bool res = false;

#ifdef _OPENMP
      omp_set_lock( &lock );
#endif
      if (head->next != tail)
      {
         Node *p = head->next;
         Node *next = p->next;
         if (item)
            *item = p->item;
         delete p;
         head->next = next;
         next->prev = head;
         res = true;
      }
#ifdef _OPENMP
      omp_unset_lock( &lock );
#endif

      return res;
   }

   bool popBack (T *item)
   {
      bool res = false;

#ifdef _OPENMP
      omp_set_lock( &lock );
#endif
      if (tail->prev != head)
      {
         Node *p = tail->prev;
         Node *prev = p->prev;
         if (item)
            *item = p->item;
         delete p;
         tail->prev = prev;
         prev->next = tail;
         res = true;
      }

#ifdef _OPENMP
      omp_unset_lock( &lock );
#endif

      return res;
   }

   bool remove (Node *q, T * item = NULL)
   {
      bool res = false;

#ifdef _OPENMP
      omp_set_lock( &lock );
#endif

      Node *p = head->next;
      while (p != tail)
      {
         if (p == q)
         {
            Node *next = p->next;
            Node *prev = p->prev;
            if (item)
               *item = p->item;
            delete p;
            prev->next = next;
            next->prev = prev;
            res = true;
            break;
         }

         p = p->next;
      }

#ifdef _OPENMP
      omp_unset_lock( &lock );
#endif
      return res;
   }

   void pushFront (const T& item, int rank = -1)
   {
      insertAfter (head, item, rank);
      //Node *p = new Node;
      //p->item = item;
      //p->rank = rank;
      //p->next = head->next;
      //p->prev = head;
      //head->next = p;
      //p->next->prev = p;
   }

   void pushBack (const T& item, int rank = -1)
   {
      insertAfter (tail->prev, item, rank);
      //Node *p = new Node;
      //p->item = item;
      //p->rank = rank;
      //p->prev = tail->prev;
      //p->next = tail;
      //tail->prev = p;
      //p->prev->next = p;
   }

   // Insert new node after node 'o'
   void insertAfter (Node *o, const T& item, int rank = -1)
   {
#ifdef _OPENMP
      omp_set_lock( &lock );
      printf("thread %d set lock item=%d rank=%d\n", omp_get_thread_num(), item, rank);
#endif

      Node *p = new Node;
      p->item = item;
      p->rank = rank;
      p->next = o->next;
      p->prev = o;
      o->next = p;
      p->next->prev = p;

#ifdef _OPENMP
      omp_unset_lock( &lock );
#endif
   }

   void insertBefore (Node *q, const T& item, int rank = -1)
   {
      insertAfter (q->prev, item, rank);
   }

   void insert (const T& item, int rank = -1)
   {
#ifdef _OPENMP
      omp_set_lock( &lock );
      printf("insert thread %d set lock item=%d rank=%d\n", omp_get_thread_num(), item, rank);
#endif
      Node *p = head->next;
      if (rank != -1)
      {
         while (p != tail and rank > p->rank)
            p = p->next;
      }

      // Default is pushFront.
      Node *node = new Node;
      node->item = item;
      node->rank = rank;
      node->prev = p->prev;
      node->next = p;
      p->prev = node;
      node->prev->next = node;
      //this->insertBefore (node, item, rank);

#ifdef _OPENMP
      omp_unset_lock( &lock );
#endif
   }

   Node* begin()
   {
      return head->next;
   }
   Node* end()
   {
      return tail;
   }
};

int main(int argc, char *argv[])
{
   int num_threads = 1;

#ifdef _OPENMP
   fprintf(stderr,"OpenMP specification c. %s %d.%d\n", TOSTRING(_OPENMP), _MY_OPENMP_MAJOR, _MY_OPENMP_MINOR);
   #ifdef _MY_OPENMP_HAS_ATOMICS
   fprintf(stderr,"OpenMP has atomics!\n");
   #endif

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

   if (argc > 1)
      n = atoi(argv[1]);

   List<int> list;

   printf("isEmpty() = %d\n", list.isEmpty());

   list.pushBack(1);

   printf("isEmpty() = %d\n", list.isEmpty());
   printf("front() = %d\n", list.front());
   printf("back() = %d\n", list.back());

   list.pushBack(2);
   printf("back() = %d\n", list.back());

   int value = -1;
   #pragma omp parallel num_threads(2) firstprivate(value)
   {
      int thread_id = 0;
#ifdef _OPENMP
      thread_id = omp_get_thread_num();
#endif
      printf("popFront(&value) = %d %d\n", list.popFront(&value), thread_id);
      printf("value = %d %d\n", value, thread_id);
      printf("front() = %d %d\n", list.front(), thread_id);
   }

   #pragma omp parallel for
   for (int i = 3; i < 7; ++i)
      list.insert(i,i);
      //list.pushFront(i,i);

   list.insert(-1,8);
   list.insert(-2,2);

   for (List<int>::Node *p = list.begin(); p != list.end(); p = p->next)
      printf("p->item = %d %d\n", p->item, p->rank);

   printf("popBack(&value) = %d\n", list.popBack(&value));
   printf("value = %d\n", value);

   for (List<int>::Node *p = list.begin(); p != list.end(); p = p->next)
      printf("p->item = %d %d\n", p->item, p->rank);

   printf("reverse\n");
   for (List<int>::Node *p = list.end()->prev; p != list.begin()->prev; p = p->prev)
      printf("p->item = %d %d\n", p->item, p->rank);

   printf("remove  %d\n", list.remove(list.begin()->next));
   for (List<int>::Node *p = list.begin(); p != list.end(); p = p->next)
      printf("p->item = %d %d\n", p->item, p->rank);

   return 0;
}
