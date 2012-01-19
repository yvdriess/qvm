#!/usr/bin/env python
##  $Id: abstfilter.py,v 1.1.1.1 2003/07/01 23:28:27 euske Exp $
##
##  abstfilter.py - A framework for cascade filters.
##


##  AbstractFeeder
##
class AbstractFeeder:

  def __init__(self, next_filter):
    self.next_filter = next_filter
    return

  def feed_next(self, s):
    self.next_filter.feed(s)
    return

  def close(self):
    self.next_filter.close()
    return


##  AbstractFilter
##
class AbstractFilter:

  def __init__(self, next_filter):
    self.next_filter = next_filter
    return

  def process(self, s):
    raise NotImplementedError
  
  def feed(self, s):
    self.feed_next(self.process(s))
    return

  def feed_next(self, s):
    self.next_filter.feed(s)
    return

  def close(self):
    self.next_filter.close()
    return


##  AbstractConsumer
##
class AbstractConsumer:

  def __init__(self):
    return

  def feed(self, s):
    raise NotImplementedError
  
  def close(self):
    return


##  FileGenerator
##
class FileGenerator(AbstractFeeder):
  def __init__(self, next_type):
    next_filter = next_type(self)
    AbstractFeeder.__init__(self, next_filter)
    self.results = []
    return
  
  def feed(self, s):
    self.results.append(s)
    return

  def close(self):
    return
  
  def pullopen(self, f):
    while 1:
      s = f.readline()
      if not s: break
      self.feed_next(s)
      if self.results:
        for s in self.results:
          yield s
        self.results = []
    for s in self.results:
      yield s
    self.results = []
    return

