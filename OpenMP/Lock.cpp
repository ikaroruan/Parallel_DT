#include "Lock.h"
#include <iostream>

Lock::Lock()
{
	omp_init_lock(&_lock);
	_is_locked = false;
	_owner = -1;
}

Lock::~Lock()
{
	omp_destroy_lock(&_lock);
}

void Lock::locked(bool value)
{
	_is_locked = value;
}

bool Lock::locked()
{
	return _is_locked;
}

bool Lock::locked_by_owner(int i)
{
	if(locked() && owner() == i)
		return true;
	
	return false;
}

void Lock::owner(int value)
{
	_owner = value;
}

int Lock::owner()
{
	return _owner;
}

bool Lock::lock()
{
	int tid = omp_get_thread_num();
	if(!locked()){
		int num = omp_test_lock(&_lock);
		if(num == 0) return false;
		locked(true);
		owner(tid);
		return true;
	}
	else{
		if(owner() == tid)
			return true;
	}

	return false;
}

/*bool Lock::lock(int* priority)
{
	int tid = omp_get_thread_num();
	if(!locked()){
		omp_set_lock(&_lock);
		locked(true);
		owner(tid);
		return true;
	}
	else{
		// Need to check priorities.
		if(priority[owner()] > priority[tid]){
			// Lock unsuccessful.
			#pragma omp critical (priority_update)
			{
			//std::cout << "Updating priority (Thread " << tid << ").\n";
			//priority[tid] = (priority[tid] + 1) % (2*omp_get_num_threads());
			//std::cout << "Done (Thread " << tid << ").\n";
			}
			return false;
		}	
		else if(priority[owner()] < priority[tid]){
			// New thread waits for locking.
			omp_set_lock(&_lock);
			locked(true);
			owner(tid);
			return true;
		}
		else{
			// Then, priority will be decided by the thread id.
			if(owner() > tid){
				// Lock unsuccessful.
				#pragma omp critical (priority_update)
				{
			//std::cout << "BUpdating priority (Thread " << tid << ").\n";
			//	priority[tid] = (priority[tid] + 1) % (2*omp_get_num_threads());
			//std::cout << "Done (Thread " << tid << ").\n";
				}
				return false;
			}
			else if(owner() < tid){
				omp_set_lock(&_lock);
				locked(true);
				owner(tid);
				return true;
			}
			else{
				// The thread already owns this lock.
				return true;
			}
		}
	}

	// Should never reach this line.
	return false;
}*/

bool Lock::lock(int* priority)
{
	bool done = false;
	bool to_return = false;
	int tid = omp_get_thread_num();
	//std::cout << "owner = " << owner() << "\n";
	//std::cout << "locked = " << (int)locked() << "\n";

	#pragma omp critical (lock_with_priority)
	{
		if(locked()){
			done = false;
		}
		else if(priority[tid] > priority[owner()])
			done = false;
		else if(priority[tid] < priority[owner()]){
			done = true;
			to_return = false;
			priority[tid] = (priority[tid] + 1) % (2*omp_get_num_threads());
		}
		else{
			if(tid == owner()){
				done = true;
				to_return = true;
			}
			else if(tid > owner()){
				done = false;
			}
			else{
				done = true;
				to_return = false;
			}
		}
	}

	if(done) return to_return;

	omp_set_lock(&_lock);
	#pragma omp critical (lock_with_priority)
	{
		locked(true);
		owner(tid);
		to_return = true;
	}
	return to_return;
}

bool Lock::unlock()
{
	if(!locked())
		return true;
	if(owner() == omp_get_thread_num()){
		locked(false);
		owner(-1);
		omp_unset_lock(&_lock);
		return true;
	}

	// In this case, the thread does not own the lock.
	return false;
}
