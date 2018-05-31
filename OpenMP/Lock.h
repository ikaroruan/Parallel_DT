#ifndef LOCK_H
#define LOCK_H

#include <omp.h>

class Lock
{
	public:
	
	Lock();
	Lock(const Lock& lk) : _lock(lk._lock), _is_locked(lk._is_locked), _owner(lk._owner) {omp_init_lock(&_lock);}
	~Lock();

	void locked(bool value);
	bool locked();
	bool locked_by_owner(int i);
	void owner(int value);
	int owner();

	bool lock();
	bool lock(int* priority);
	bool unlock();

	private:

	omp_lock_t _lock;
	bool _is_locked;
	int _owner;
};

#endif // LOCK_H
