#include "LPdefine.h"

int HashValue(hashtab *tab, char *key, unsigned scal)
{
  unsigned val;
  char *s;

  val = 0;
  for (s = key; *s != '\0'; s++)
    val = scal * val + (*s);

  return (val % tab->lssz);
} /* HashValue */

int InsertHash(hashtab *tab, char *key, int addr, unsigned scal)
{
  hashlst *ptr;
  int k, len;

  k = HashValue(tab, key, scal);

  for (ptr = tab->list[k]; ptr; ptr = ptr->next)
    if (!CompString(ptr->key, key))
      break;

  if (!ptr)
  {
    ptr = (hashlst *)cAlloc(sizeof(hashlst),
                            "ptr, InsertHash");

    len = strlen(key);
    ptr->key = cAlloc(len + 1,
                      "for ptr->key in InsertHash");

    strcpy(ptr->key, key);
    ptr->addr = addr;

    ptr->next = tab->list[k];
    tab->list[k] = ptr;
    return false;
  }

  return true;
} /* InsertHash */

int FoundHash(hashtab *tab, char *key, unsigned scal)
{
  int k, addr;
  hashlst *ptr;

  addr = -1;
  k = HashValue(tab, key, scal);

  for (ptr = tab->list[k]; ptr; ptr = ptr->next)
  {
    if (!CompString(ptr->key, key))
    {
      addr = ptr->addr;
      break;
    }
  }

  return addr;
} /* FoundHash */
