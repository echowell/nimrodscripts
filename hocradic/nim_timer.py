#!/usr/bin/env python3
import time
import functools

class timeObject:
    def __init__(self,name):
        self.__start=None
        self.__elapsed=0.0
        self.__name=name
      
    def start(self):
        if self.__start is not None:
            print(f"{self.__name} timer was not reset")
            raise ValueError
        self.__start=time.time()
        return None
      
    def stop(self):
        end=time.time()
        if self.__start is None:
            print(f"{self.__name} timer was never started")
            raise ValueError
        self.__elapsed+=(end-self.__start)
        self.__start=None
        return None
      
    def get_elpased(self):
        return self.__elapsed
      
    def print_time(self):
        print(f"{self.__name}: {self.__elapsed}s")
        return None
    
class nimTimer:
    ''' A timer class for profiling python scripts '''
    def __init__(self):
        self.__ids={}
    
    def start(self,id):
        if id in self.__ids:
            self.__ids[id].start()
        else:
            self.__ids[id]=timeObject(id)
            self.__ids[id].start()
        return None
      
    def stop(self,id):
        if id not in self.__ids:
            print(f"{id} not found in nimTime")
            raise Keyerror
        self.__ids[id].stop()
        return None

    def get_elapsed(self,id):
        if id not in self.__ids:
            print(f"{id} not found in nimTime")
            raise Keyerror
        return self.__ids[id].get_elapsed()
      
    def print_times(self):
        for val in self.__ids.values():
            val.print_time()
        return None

timer=nimTimer()

def timer_func(func):
    '''function wrapper to use nimTimer as a decorator
       use: @nim_timer.timer_func before function def
    '''
    @functools.wraps(func)
    def wrapper(*args,**kwargs):
        timer.start(func.__name__)
        result=func(*args,**kwargs)
        timer.stop(func.__name__)
        return result
    return wrapper
