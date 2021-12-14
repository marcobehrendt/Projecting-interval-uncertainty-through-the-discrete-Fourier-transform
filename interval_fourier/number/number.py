"""
    #---------------------------------#
    Created on Tue Dec 22 10:58:17 2020
    @author: Marco De Angelis
    University of Liverpool
    github.com/marcodeangelis
    GNU LGPL v3.0
    #---------------------------------#

    Number module for the *Fourier Transform* library.
"""
import numpy
from itertools import repeat
import matplotlib.pyplot as pyplot # used in IntervalVector only

def intervalDataTypes(): # data types (or classes) in this module
    return ('Interval','ComplexInterval','IntervalVector') # Data types sharing same public methods, like lo(), hi(), mid(), etc...
def numpyDataTypes():
    return ('int8','int16','int32','int64','float16','float32','float64','float_','complex128','complex_') # not exhaustive
def numberDataTypes():
    return ('int','float','complex','int8','int16','int32','int64','float16','float32','float64','complex128') # not exhaustive
def complexNumberTypes():
    return(['complex128','complex_','complex'])

# The code in this file should comply to PEP-8: https://realpython.com/python-pep8/

machine_epsilon = 7./3 - 4./3 - 1

class Interval():                  # simple interval class
    """
    Created on Tue Dec 26 11:59:25 2017
    @author: Marco De Angelis
    University of Liverpool
    github.com/marcodeangelis
    GNU LGPL v3.0
    """
    def __repr__(self): # return
        return "【%g, %g】"%(self.__lo,self.__hi)  # https://www.compart.com/en/unicode/U+3011 #lenticular brackets

    def __str__(self): # print
        return "【%g, %g】"%(self.__lo,self.__hi)  # https://www.compart.com/en/unicode/U+3011 #lenticular brackets

    def __init__(self,*args):
        if (args is None) | (len(args)==0):
            self.__lo, self.__hi = -1, 1
        if len(args)==1:
            if args[0].__class__.__name__ in numberDataTypes():
                self.__lo, self.__hi = args[0], args[0]
            elif args[0].__class__.__name__ == 'Interval':
                self = args[0]
            else:
                self.__lo, self.__hi = args[0][0], args[0][1]
        if len(args)==2:
            self.__lo, self.__hi = args[0], args[1]
    ## Class methods start here ##
    def __hash__(self):       #Makes interval class hashable     # https://docs.python.org/3/reference/datamodel.html#object.__hash__
        return hash((self.lo(),self.hi()))
    def lo(self):
        return self.__lo
    def hi(self):
        return self.__hi
    def mid(self):
        return (self.__lo + self.__hi)/2
    def rad(self):
        return (self.__hi - self.__lo)/2
    def width(self):
        return self.__hi - self.__lo
    def inf(self): #  != lo. support for outword directed rounding
        return self.__lo - machine_epsilon
    def sup(self): #  != hi. support for outword directed rounding
        return self.__hi + machine_epsilon
    def stradzero(self): # iszeroin
        if (self.__lo <= 0) & (self.__hi >= 0): return True 
        else: return False
    def contains(self,other): # True if self contains other
        if other.__class__.__name__ not in intervalDataTypes():
            other = Interval(other,other)
        return (self.inf() <= other.inf()) & (self.sup() >= other.sup())
    def encloses(self,other): # True if self encloses other (strictly on both sides)
        return self.lo() < other.lo() and other.hi() < self.hi()
    def inside(self,other): # True if other contains self
        if other.__class__.__name__ not in intervalDataTypes(): # other is a scalar
            other = Interval(other,other)
        return (self.inf() >= other.inf()) & (self.sup() <= other.sup())
    def intersect(self,other): # True if self intersects other and viceversa
        return not(self < other or other < self)
    def union(self,other): # Some say union can only be done between intersecting intervals.
        return Interval(min(self.lo(),other.lo()),max(self.hi(),other.hi()))
    def intersection(self,other):
        if self.intersect(other):
            return Interval(max(self.lo(),other.lo()), min(self.hi(),other.hi()))
        else:
            return None
    def slider(self,p):
        if p.__class__.__name__ in ['list','tuple']:
            return [self.__lo + pi * self.width() for pi in p]
        else:
            return self.__lo + p * self.width()
    def linspace(self,N=30):
        return list(numpy.linspace(self.lo(),self.hi(),num=N))
    ## Class methods end here ##

    #-------------------------------------#
    # Override arithmetic operators START #
    #-------------------------------------#
    # unary operators #
    def __neg__(self):
        return Interval(-self.__hi, -self.__lo)
    def __pos__(self):
        return self
    # binary operators #
    def __add__(self,other):
        otherType = other.__class__.__name__
        if otherType in numberDataTypes(): # add support for numpy
            if otherType in complexNumberTypes():
                addL = self.__lo + other
                addH = self.__hi + other
                return ComplexInterval(addL,addH)
            else:
                addL = self.__lo + other
                addH = self.__hi + other
                return Interval(addL,addH)
        elif otherType in intervalDataTypes():
            if otherType == 'ComplexInterval':
                addL = self.__lo + other.lo()
                addH = self.__hi + other.hi()
                return ComplexInterval(addL,addH)
            else:
                addL = self.__lo + other.__lo
                addH = self.__hi + other.__hi
                return Interval(addL,addH)
        else:
            raise TypeError('Addition only allowed between numbers.')
    def __radd__(self, left):
        leftType = left.__class__.__name__
        if leftType in numberDataTypes():
            if leftType in complexNumberTypes():
                addL = self.__lo + left
                addH = self.__hi + left
                return ComplexInterval(addL,addH)
            else:
                addL = left + self.__lo
                addH = left + self.__hi
                return self.__add__(left)
        else:
            raise TypeError('Addition only allowed between numbers.')
    def __sub__(self, other):
        otherType = other.__class__.__name__
        if otherType in numberDataTypes():
            if otherType in complexNumberTypes(): # needs extending to numpy complex types too.
                subL = self.__lo - other
                subH = self.__hi - other
                return ComplexInterval(subL,subH)
            else:
                subL = self.__lo - other
                subH = self.__hi - other
                return Interval(subL,subH)
        elif otherType in intervalDataTypes():
            if otherType == 'ComplexInterval':
                subL = self.__lo - other.lo()
                subH = self.__hi - other.hi()
                return ComplexInterval(subL,subH)
            else:
                subL = self.__lo - other.__hi
                subH = self.__hi - other.__lo
                return Interval(subL,subH)
    def __rsub__(self, left):
        leftType = left.__class__.__name__
        if leftType in numberDataTypes():
            if leftType in complexNumberTypes():
                subL = self.__lo - left
                subH = self.__hi - left
                return ComplexInterval(subL,subH)
            else:
                subL = left - self.__hi
                subH = left - self.__lo
                return Interval(subL,subH)
        else:
            raise TypeError('Subtraction only allowed between numbers.')
    def __mul__(self,other):
        otherType = other.__class__.__name__
        if otherType in numberDataTypes():
            if otherType in complexNumberTypes():
                Real = self*other.real - 0*other.imag
                Imag = self*other.imag + 0*other.real
                mulL = Real.lo()+1j*Imag.lo()
                mulH = Real.hi()+1j*Imag.hi()
                return ComplexInterval(mulL,mulH)
            else:
                if other>0:
                    mulL = self.__lo * other
                    mulH = self.__hi * other
                else:
                    mulL = self.__hi * other
                    mulH = self.__lo * other
        elif otherType in intervalDataTypes():
            if otherType == 'ComplexInterval':
                return other*self # * will be done in ComplexInterval
            else:
                if (self.__lo>=0) & (other.__lo>=0): # A+ B+
                    mulL = self.__lo * other.__lo
                    mulH = self.__hi * other.__hi
                elif (self.__lo>=0) & ((other.__lo<0) & (other.__hi>0)): # A+ B0
                    mulL = self.__hi * other.__lo
                    mulH = self.__hi * other.__hi
                elif (self.__lo>=0) & (other.__hi<=0): # A+ B-
                    mulL = self.__hi * other.__lo
                    mulH = self.__lo * other.__hi
                elif ((self.__lo<0) & (self.__hi>0)) & (other.__lo>=0): # A0 B+
                    mulL = self.__lo * other.__hi
                    mulH = self.__hi * other.__hi
                elif ((self.__lo<0) & (self.__hi>0)) & ((other.__lo<0) & (other.__hi>0)): # A0 B0
                    mulL1 = self.__lo * other.__hi
                    mulL2 = self.__hi * other.__lo
                    mulL = min(mulL1,mulL2)
                    mulH1 = self.__lo * other.__lo
                    mulH2 = self.__hi * other.__hi
                    mulH = max(mulH1,mulH2)
                elif ((self.__lo<0) & (self.__hi>0)) & (other.__hi<=0): # A0 B-
                    mulL = self.__hi * other.__lo
                    mulH = self.__lo * other.__lo
                elif (self.__hi<=0) & (other.__lo>=0): # A- B+
                    mulL = self.__lo * other.__hi
                    mulH = self.__hi * other.__lo
                elif (self.__hi<=0) & ((other.__lo<0) & (other.__hi>0)): # A- B0
                    mulL = self.__lo * other.__hi
                    mulH = self.__lo * other.__lo
                elif (self.__hi<=0) & (other.__hi<=0): # A- B-
                    mulL = self.__hi * other.__hi
                    mulH = self.__lo * other.__lo
        return Interval(mulL,mulH)
    def __rmul__(self, left):
        if left.__class__.__name__ in numberDataTypes():
            return self.__mul__(left)
        else:
            raise TypeError('Multiplication only allowed between numbers.')
    def __truediv__(self,other):
        if other.__class__.__name__ in numberDataTypes():
            if other>0:
                divL = self.__lo / other
                divH = self.__hi / other
            elif other<0:
                divL = self.__hi / other
                divH = self.__lo / other
        elif other.__class__.__name__ in intervalDataTypes():
            if other.stradzero():
                raise Warning("Division by interval containing zero") # TODO: extended interval arithmetic
            if (self.__lo>=0) & (other.__lo>0):
                divL = self.__lo/other.__hi
                divH = self.__hi/other.__lo
            elif ((self.__lo<0) & (self.__hi>0)) & (other.__lo>0):
                divL = self.__lo/other.__lo
                divH = self.__hi/other.__lo
            elif (self.__hi<=0) & (other.__lo>0):
                divL = self.__lo/other.__lo
                divH = self.__hi/other.__hi
            elif (self.__lo>=0) & (other.__hi<0):
                divL = self.__hi/other.__hi
                divH = self.__lo/other.__lo
            elif ((self.__lo<0) & (self.__hi>0)) & (other.__hi<0):
                divL = self.__hi/other.__hi
                divH = self.__lo/other.__hi
            elif (self.__hi<=0) & (other.__hi<0):
                divL = self.__hi/other.__lo
                divH = self.__lo/other.__hi
        return Interval(divL,divH)
    def __rtruediv__(self, left):
        if left.__class__.__name__ in numberDataTypes():
            if left>0:
                if (self.__lo>0):
                    divL = left / self.__hi
                    divH = left / self.__lo
                elif (self.__hi<0):
                    divL = left / self.__hi
                    divH = left / self.__lo
                else:
                    raise ZeroDivisionError('Division is allowed for intervals not containing the zero') # this should not return an error, but rather an unbounded interval
            elif left<0:
                if (self.__lo>0):
                    divL = left / self.__lo
                    divH = left / self.__hi
                elif (self.__hi<0):
                    divL = left / self.__lo
                    divH = left / self.__hi
                else:
                    raise ZeroDivisionError('Division is allowed for intervals not containing the zero') # this should not return an error, but rather an unbounded interval
            return Interval(divL,divH)
        else:
            raise TypeError('Division only allowed between numbers.')
    def __pow__(self,other):
        otherType = other.__class__.__name__
        if otherType in intervalDataTypes():
            raise TypeError('Power elevation by an interval not needed for this library.')
        elif otherType in numberDataTypes():
            if (other%2==0) | (other%2==1):
                other = int(other)
            if otherType == "int":
                if other > 0:
                    if other%2 == 0: # even power
                        if self.__lo >= 0:
                            powL = self.__lo ** other
                            powH = self.__hi ** other
                        elif self.__hi < 0:
                            powL = self.__hi ** other
                            powH = self.__lo ** other
                        else: # interval contains zero
                            H = max(-self.__lo,self.__hi)
                            powL = 0
                            powH = H ** other
                    elif other%2 == 1: # odd power
                        powL = self.__lo ** other
                        powH = self.__hi ** other
                elif other < 0:
                    if other%2 == 0: # even power
                        if self.__lo >= 0:
                            powL = self.__hi ** other
                            powH = self.__lo ** other
                        elif self.__hi < 0:
                            powL = self.__lo ** other
                            powH = self.__hi ** other
                        else: # interval contains zero
                            print("Error. \nThe interval contains zero, so negative powers should return \u00B1 Infinity")
                    elif other%2 == 1: # odd power
                        if self.__lo != 0:
                            powL = self.__hi ** other
                            powH = self.__lo ** other
                        else: # interval contains zero
                            print("Error. \nThe interval contains zero, so negative powers should return \u00B1 Infinity")
            elif otherType == "float":
                    if self.__lo >= 0:
                        if other > 0:
                            powL = self.__lo ** other
                            powH = self.__hi ** other
                        elif other < 0:
                            powL = self.__hi ** other
                            powH = self.__lo ** other
                        elif other == 0:
                            powL = 1
                            powH = 1
        return Interval(powL,powH)
    def __rpow__(self,left):
        raise TypeError('Interval exponents are not needed for this library.')
    def __lt__(self, other):
        if other.__class__.__name__ in intervalDataTypes():
            return self.sup() < other.inf()
        elif other.__class__.__name__ in numberDataTypes():
            return self.sup() < other
    def __rlt__(self,left):
        if left.__class__.__name__ in intervalDataTypes():
            return left.sup() < self.inf()
        elif left.__class__.__name__ in numberDataTypes():
            return left < self.inf()
    def __gt__(self, other):
        if other.__class__.__name__ in intervalDataTypes():
            return self.inf() > other.sup()
        elif other.__class__.__name__ in numberDataTypes():
            return self.inf() > other
    def __rgt__(self, left):
        if left.__class__.__name__ in intervalDataTypes():
            return left.inf() > self.sup()
        elif left.__class__.__name__ in numberDataTypes():
            return left > self.sup()
    def __le__(self, other):
        if other.__class__.__name__ in intervalDataTypes():
            return self.sup() <= other.inf()
        elif other.__class__.__name__ in numberDataTypes():
            return self.sup() <= other
    def __rle__(self,left):
        if left.__class__.__name__ in intervalDataTypes():
            return left.sup() <= self.inf()
        elif left.__class__.__name__ in numberDataTypes():
            return left <= self.inf()
    def __ge__(self, other):
        if other.__class__.__name__ in intervalDataTypes():
            return self.inf() >= other.sup()
        elif other.__class__.__name__ in numberDataTypes():
            return self.inf() >= other
    def __rge__(self, left):
        if left.__class__.__name__ in intervalDataTypes():
            return left.inf() >= self.sup()
        elif left.__class__.__name__ in numberDataTypes():
            return left >= self.sup()
    def __eq__(self, other):
        if other.__class__.__name__ in intervalDataTypes():
            return hash(self)==hash(other)
        else:
            return False
    def __ne__(self,other):
        return not(self == other)
    #------------------------------------------------------------------------------------------------------
    # Override arithmetic operators END
    #------------------------------------------------------------------------------------------------------


class ComplexInterval():            # simple complex interval class
    """
    Created on Mon Jul 13 17:02:25 2020
    @author: Marco De Angelis
    University of Liverpool
    github.com/marcodeangelis
    GNU LGPL v3.0

    Within the *gappy Fourier transform* library this class is only used for addition and multiplication.
    Support for subintervalization was also removed to make the code slimmer.
    """
    def __repr__(self): # return
        return "【%g%+gi, %g%+gi】"%(self.__lo.real,self.__lo.imag,self.__hi.real,self.__hi.imag)

    def __str__(self): # print
        return "【%g%+gi, %g%+gi】"%(self.__lo.real,self.__lo.imag,self.__hi.real,self.__hi.imag)

    def __init__(self,*args):
        if (args is None) | (len(args)==0):
            lo, hi = -1-1j, 1+1j
        if len(args)==1:
            lo, hi = args[0], args[0]
        if len(args)==2:
            lo, hi = args[0], args[1]
        self.__lo = lo
        self.__hi = hi
    # ------------------------------------ #
    ## ---- Class methods start here ---- ##
    def lo(self):
        return self.__lo
    def hi(self):
        return self.__hi
    def mid(self):
        return (self.__lo + self.__hi)/2
    def rad(self):
        return (self.__hi - self.__lo)/2
    def width(self):
        return self.__hi - self.__lo
    def stradzero(self):
        iszeroin = [False,False]
        if (self.__lo.real <= 0) & (self.__hi.real >= 0):
            iszeroin[0] = True
        if (self.__lo.imag <= 0) & (self.__hi.imag >= 0):
            iszeroin[1] = True
        return iszeroin
    def slider(self,p):
        return self.__lo + p * self.width()
    def real(self):
        return Interval(self.__lo.real, self.__hi.real)
    def imag(self):
        return Interval(self.__lo.imag, self.__hi.imag)
    def conjugate(self):
        return ComplexInterval(self.__lo.conjugate(),self.__hi.conjugate())
    def absolute(self):
        return (self.real()**2 + self.imag()**2)**0.5
    def __abs__(self):
        return self.absolute()
    ## ------ Class methods end here ------ ##
    #----------------------------------------#
    # Override arithmetical operations START #
    #----------------------------------------#
    def __add__(self,other):
        otherType = other.__class__.__name__
        if otherType in numberDataTypes(): # add support for numpy
            addL = self.__lo + other
            addH = self.__hi + other
            return ComplexInterval(addL,addH)
        elif otherType in intervalDataTypes():
            addL = self.__lo + other.lo()
            addH = self.__hi + other.hi()
            return ComplexInterval(addL,addH)
        else:
            raise TypeError('Addition only allowed between numbers.')
    def __radd__(self, left):
        leftType = left.__class__.__name__
        if leftType in numberDataTypes():
            return self.__add__(left)
        else:
            raise TypeError('Addition only allowed between numbers.')
    def __sub__(self, other):
        otherType = other.__class__.__name__
        if otherType in numberDataTypes():
            subL = self.__lo - other
            subH = self.__hi - other
            return ComplexInterval(subL,subH)
        elif otherType in intervalDataTypes():
            subL = self.__lo - other.hi()
            subH = self.__hi - other.lo()
            return ComplexInterval(subL,subH)
    def __rsub__(self, left):
        leftType = left.__class__.__name__
        if leftType in numberDataTypes():
            subL = left - self.__hi
            subH = left - self.__lo
            return ComplexInterval(subL,subH)
        else:
            raise TypeError('Subtraction only allowed between numbers.')
    def __mul__(self,other):
        otherType = other.__class__.__name__
        if otherType in numberDataTypes():
            Real = self.real()*other.real - self.imag()*other.imag
            Imag = self.real()*other.imag + self.imag()*other.real
            mulL = Real.lo()+1j*Imag.lo()
            mulH = Real.hi()+1j*Imag.hi()
            return ComplexInterval(mulL,mulH)
        elif otherType in intervalDataTypes():
            if otherType == 'ComplexInterval':
                Real = self.real()*other.real() - self.imag()*other.imag()
                Imag = self.real()*other.imag() + self.imag()*other.real()
                mulL = Real.lo()+1j*Imag.lo()
                mulH = Real.hi()+1j*Imag.hi()
                return ComplexInterval(mulL,mulH)
            else:
                Real = self.real()*other # - self.imag()*other.imag()
                Imag = self.imag()*other #self.real()*other.imag() + self.imag()*other.real()
                mulL = Real.lo()+1j*Imag.lo()
                mulH = Real.hi()+1j*Imag.hi()
                return ComplexInterval(mulL,mulH)
    def __rmul__(self, left):
        leftType = left.__class__.__name__
        if leftType in numberDataTypes():
            return self.__mul__(left)
        else:
            raise TypeError('Multiplication only allowed between numbers.')
    def __truediv__(self,other):
        otherType = other.__class__.__name__
        if otherType in numberDataTypes():
            a,b,c,d = self.real(), self.imag(), other.real, other.imag
            Real = (a*c + b*d)/(c**2 + d**2)
            Imag = (b*c - a*d)/(c**2 + d**2)
            divL = Real.lo()+1j*Imag.lo()
            divH = Real.hi()+1j*Imag.hi()
            return ComplexInterval(divL,divH)
        elif otherType in intervalDataTypes():
            if otherType == 'ComplexInterval':
                a,b = self.real(), self.imag()#, other.real(), other.imag()
                a,b,c,d = self.real(), self.imag(), other.real(), other.imag()
                Real = (a*c + b*d)/(c**2 + d**2)
                Imag = (b*c - a*d)/(c**2 + d**2)
                divL = Real.lo()+1j*Imag.lo()
                divH = Real.hi()+1j*Imag.hi()
                return ComplexInterval(divL,divH)
            else:
                a,b,c = self.real(), self.imag(), other
                Real = a/c
                Imag = b/c
                divL = Real.lo()+1j*Imag.lo()
                divH = Real.hi()+1j*Imag.hi()
                return ComplexInterval(divL,divH)
    def __rtruediv__(self, left):
        leftType = left.__class__.__name__
        if leftType in numberDataTypes():
            if leftType in complexNumberTypes():
                a,b,c,d = left.real, left.imag, self.real(), self.imag()
                Real = (a*c + b*d)/(c**2 + d**2)
                Imag = (b*c - a*d)/(c**2 + d**2)
                divL = Real.lo()+1j*Imag.lo()
                divH = Real.hi()+1j*Imag.hi()
                return ComplexInterval(divL,divH)
            else:
                a,c,d = left, self.real(), self.imag()
                Real = (a*c)/(c**2 + d**2)
                Imag = -(a*d)/(c**2 + d**2)
                divL = Real.lo()+1j*Imag.lo()
                divH = Real.hi()+1j*Imag.hi()
                return ComplexInterval(divL,divH)
        else:
            raise TypeError('Division only allowed between numbers.')
    #-------------------------------------#
    # Override arithmetical operations END
    #-------------------------------------#

class IntervalVector(): # wrapper class of the scalar class interval with some plotting facilities
    """
    Created on Mon Jul 24 17:09:51 2020
    @author: Marco De Angelis
    University of Liverpool
    github.com/marcodeangelis
    GNU LGPL v3.0
    """
    def __repr__(self): # return
        if len(self)>10:
            a = [str(i) for i in self]
            s = '\n'.join(a[:5]+['...']+a[-5:-1])
        else:
            s = '\n'.join([str(i) for i in self])
        return s
    def __str__(self): # print
        return self.__repr__()
    def __len__(self):
        return len([i for i in self])
    def __init__(self,*args,notation='infsup',axis=1,name=''):
        self.name = name # initialise this with an empty string
        if len(args)==0:  # what should an empty IntervalArray(object) look like?
            self.__lo = [-1]
            self.__hi = [1]
        elif len(args)==1:   # this must be a list, tuple (array?) of intervals
            assert args[0].__class__.__name__ in ['list', 'tuple','IntervalVector'], 'single input must be list or a tuple of intervals.'
            if args[0][0].__class__.__name__ in ['Interval','ComplexInterval']:
                self.__lo = [x.lo() for x in args[0]]
                self.__hi = [x.hi() for x in args[0]]
            elif args[0][0].__class__.__name__ in ['list','tuple']:
                if axis == 0:
                    assert len(args[0]) == 2, 'a list or tuple is needed.'
                    self.__lo, self.__hi = args[0][0], args[0][1]
                elif axis == 1:
                    self.__lo = list([x for x in zip(*args[0])][0])
                    self.__hi = list([x for x in zip(*args[0])][1])
        elif len(args)==2:
            if args[0].__class__.__name__ == 'list':
                self.__lo, self.__hi = args[0], args[1]
            else:
                self.__lo = list()
                self.__hi = list()
                for a in args:
                    if a.__class__.__name__ == 'tuple':
                        self.__lo.append(a[0]) 
                        self.__hi.append(a[1])
                    elif a.__class__.__name__ == 'Interval':
                        self.__lo.append(a.lo()) 
                        self.__hi.append(a.hi())
                    else:
                        raise TypeError('multiple arguments must be a tuple or an interval.')
    def __iter__(self): # make class iterable
        for l,u in zip(self.__lo,self.__hi):
            yield Interval(l,u)
    def __getitem__(self,index): # make class subscrictable
        if index.__class__.__name__ in ['list','tuple']:
            if len(index)>0:
                return IntervalVector([Interval(self.__lo[i],self.__hi[i]) for i in index])
            else:
                return IntervalVector([]) # todo: create empty dataset
        else:
            return Interval(self.__lo[index],self.__hi[index])
    def inf(self):
        return self.__lo
    def lo(self):
        return self.__lo
    def sup(self):
        return self.__hi
    def hi(self):
        return self.__hi
    def tolist(self):
        return [Interval(l,h) for l,h in zip(self.__lo,self.__hi)]
    def toarray(self, order='F'):
        if order=='F':
            return numpy.array([self.__lo, self.__hi])
        elif order=='C':
            return numpy.array([self.__lo, self.__hi]).T
    def slider(self,p=0.5):  
        if p.__class__.__name__ in ['list','tuple']:
            assert len(self)==len(p), f'p must be of length {len(self)}'
            return [si.slider(pi) for si,pi in zip(self,p)]
        else:
            return [si.slider(p) for si in self] # p = list(repeat(.5, times=len(self)))
    def rand(self,N=1):
        n = len(self)
        r = numpy.random.random_sample(size=(N,n))
        lo_arr = numpy.array(N*[self.lo()])
        hi_arr = numpy.array(N*[self.hi()])
        if N>1:
            in_arr = lo_arr + r * (hi_arr-lo_arr)
        elif N==1:
            in_arr = (lo_arr + r * (hi_arr-lo_arr))[0]
        return in_arr
    # Magic methods
    def __add__(self,other):
        return IntervalVector([a+b for a,b in zip(self,other)])
    def __sub__(self,other):
        return IntervalVector([a-b for a,b in zip(self,other)])
    def __mul__(self,other):
        return IntervalVector([a*b for a,b in zip(self,other)])
    def __truediv__(self,other):
        return IntervalVector([a/b for a,b in zip(self,other)])
    def plot(self,marker='_',size=20,xlabel='x',ylabel='y',title='',save=None,ax=None,label=None,alpha=0.3):
        N = len(self)
        if ax is None:
            fig = pyplot.figure(figsize=(18,6))
            ax = fig.subplots()
            ax.grid()
        x = list(range(0,N))
        ax.plot(x,self.lo())
        ax.plot(x,self.hi())
        ax.fill_between(x=x, y1=self.hi(), y2=self.lo(), alpha=alpha,label=label)
        for i in range(0,N):
            ax.scatter([i,i],[self.lo()[i],self.hi()[i]],s=size,marker=marker)
        ax.set_xlabel(xlabel,fontsize=20)
        ax.set_ylabel(ylabel,fontsize=20)
        ax.tick_params(direction='out', length=6, width=2, labelsize=14) #https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.axes.Axes.tick_params.html
        ax.set_title(title,fontsize=20)
        if save is not None:
            ax.savefig(save)
        return None