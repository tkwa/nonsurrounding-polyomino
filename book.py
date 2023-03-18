#math functions to use on a whim
#Add comments, future-me; you're not self-omniscient
#http://www.sympy.org/en/features.html for future maybe
import math
import random
def gcd(*args):
    """
    Greatest common divisor of positive integers,
    or of the list if one argument is given.
    """
    if(len(args)==0):
        return 1
    elif(len(args)==1):
        l=args[0]
        if(len(l)==2):
            return gcd(l[0],l[1])
        return gcd(l[0],gcd(l[1:]))
    elif(len(args)==2):
        if args[0] == 0:
            return args[1]
        return gcd(args[1] % args[0], args[0])
    else:
        return gcd(args[0],gcd(args[1:]))      
def lcm(*args):
    """
    Least common multiple of positive integers,
    or of the list if one argument is given.
    """
    if(len(args)==0):
        return 1
    elif(len(args)==1):
        l=args[0]
        if(len(l)==0):
            return 1
        if(len(l)==2):
            return lcm(l[0],l[1])
        return lcm(l[0],lcm(l[1:]))
    elif(len(args)==2):
        return args[0]*args[1]//gcd(args[0],args[1])
    else:
        return lcm(args[0],lcm(args[1:]))
def sieve(n):
    """
    Taken from
    https://stackoverflow.com/questions/16004407/
    a-fast-prime-number-sieve-in-python
    """
    sz = n//2
    s = [1]*sz
    limit = int(n**0.5)
    for i in range(1,limit):
        if s[i]:
            val = 2*i+1
            tmp = ((sz-1) - i)//val 
            s[i+val::val] = [0]*tmp
    return [2] + [i*2+1 for i, v in enumerate(s) if v and i>0]
primes=sieve(10000)
def modpow(a,e,n):
    """Returns a^e (mod n).
    More efficient for large values than directly computing."""
    if(n<1):
        return -1
    b=bin(e)[2:]
    prod=1
    current=a
    for i in range(len(b)):
        if(b[-i-1]=='1'):
            prod *= current
            prod %= n
        current *= current
        current %= n
    return prod
initialCertainty=1 #Value to get approximation for Bayesian probability a
#number given by Miller-Rabin is prime
for p in primes:
        initialCertainty*=1+1/p
def ic_hack():
    return initialCertainty
def v(n,d):
    """Maximum e such that d^e divides n."""
    if n<1 or d<2:
        return -1
    count=0
    while n%d==0:
        n//=d
        count+=1
    return count
def millerRabin(n,prob=-1,checked=-1):
    """Probabilistic primality test. Defaults to 99.999% Bayesian expectation
    with priors given by PNT but can be set higher with prob.
    You may include an optional array of primes which are known
    not to divide the integer (defaults to the primes <10000
    for the sake of the standard book.prime)."""
    if prob==-1:
        prob=0.99999
    initialCertainty=ic_hack()
    if checked != -1:
        initialCertainty=1
        for p in checked:
            initialCertainty*=1+1/p
    desired=1/(1-prob)
    certainty=1/math.log(n)
    certainty/=1-certainty
    certainty*=initialCertainty
    
    s=v(n-1,2)
    d=(n-1)//2**s
    count=0
    while certainty<desired or count<10:
        a=random.randint(1,n-1)
        k=modpow(a,d,n)
        if k != 1:
            r=0
            while k%n!=n-1 and r<s:
                k=(k**2)%n
                r+=1
            if r==s:
                return False
        certainty*=4
        count+=1
    return True
def prime(n,prob=-1):
    """Deterministic up to 100,000,000,
    uses Miller-Rabin afterwards to Bayesian probability of at least prob.
    Default value set to 99.999%.
    """
    if(n<2):
        return False
    if(n==2):
        return True
    for p in primes:
        if (n%p==0):
            return False
        elif (p*p>n):
            return True
    return millerRabin(n,prob)
def dprime(n):
    """Deterministic (and slow) primality test. May run for a long time."""
    if(n<2):
        return False
    if(n==2):
        return True
    for p in primes:
        if (n%p==0):
            return False
        elif (p*p>n):
            return True
    i=primes[len(primes)-1]+2 
    while(i*i<=n):
        if(n%i==0):
            return False
        i+=2
    return True
def p(n):
    "nth prime"
    ps=primes
    if(n<=len(primes)):
        return primes[n-1]
    primecount=len(primes)
    t=primes[-1]+2
    while(primecount<n):
        if(prime(t)):
            primecount+=1
        t+=2
    return t-2
def addprimepower(m,pp,i=0):
    """Takes in a factored m and an array [prime,exp]
    and makes m's factorization multiply by prime**exp.
    i is a minimum index of the prime in m, if known"""
    n=m[:]
    if pp[1]>0:
        j=i
        while j<len(n):
            if n[j][0]==pp[0]:
                n[j][1]+=pp[1]
                return n
            elif n[j][0]>pp[0]:
                return n[:j]+[pp]+n[j:]
            j+=1
        return n+[pp]
    else:
        return n
def divmult(*args):
    """Takes in a number of factorizations and returns
    a factorization of their product."""
    d={}
    for l in args:
        for e in l:
            if e[0] not in d:
                d[e[0]]=e[1]
            else:
                d[e[0]]+=e[1]
    return [[p,d[p]] for p in d]
def divdiv(a,b):
    """Takes in the factorizations of two integers and returns the
    factorization of their quotient. Fails if b does not divide a."""
    d={e[0]:e[1] for e in a}
    for e in b:
        if e[0] in d and e[1]<=d[e[0]]:
            d[e[0]]-=e[1]
            if d[e[0]]==0:
                del d[e[0]]
        else:
            return -1
    return [[p,d[p]] for p in d]
def dsieve(n,s=[]):
    """Returns an array of all the integers up to n
    in factored form. More efficient than computing directly.
    Takes in an array of the primes up to n if already available."""
    if s==[]:
        s=sieve(n)
    out=[]
    for i in range(n):
        out.append([])
    for p in s:
        k=p
        while k<n:
            out[k].append([p,v(k,p)])
            k+=p
    return out
def div(n,intermediate=0):
    """
    Factorization of n, given in an array of arrays of the form [prime,exponent]
    in increasing order of prime size
    E.g. div(24)=[[2,3],[3,1]]
    Provide an array of numbers which multiply to give n
    as the second argument to speed up factorization.
    """
    if n==0:
        return []
    if intermediate!=0:
        result=[]
        prod=1
        for p in intermediate:
            result=divmult(result,div(p))
            prod*=p
        if prod==n:
            return result
        else: #if the factorization was wrong
            return div(n)
    out=[]
    if(n>1):
        for p in primes:
            if(n%p==0):
                out.append([p,0])
                while(n%p==0):
                    n=n//p
                    out[-1][1]+=1
        j=primes[-1]+2
        if(prime(n)):
            out+=[[n,1]]
            return out
        while(n>1):
            if(n%j==0):
                out.append([j,0])
                while(n%j==0):
                    n=n//j
                    out[-1][1]+=1
            j+=2
            if(j%3==0 or j%5==0):
                j+=2
        return out
    else:
        return []
    ###
    ### Many of the following functions have an optional final argument for
    ### the divisors of the input; if the program reuses a factorization,
    ### this can be very time-saving.
    ###
def undiv(d):
    p=1
    for e in d:
        p*=e[0]**e[1]
    return p
def unsorted_factors(n,d=0):
    """
    Returns an unsorted array of the factors of n.
    Do not use for things like sigma(n),
    which are significantly faster to call directly.
    """
    if(d==0):
        d=div(n)
    if(n==1):
        return [1]
    f=[]
    track=[0]*len(d)
    for i in range(0,len(d)):
        track[i]=d[i][1]
    track[len(d)-1]+=1
    while(track!=[0]*len(d)):
        ind=len(d)-1
        done=False
        while(ind>=0 and done==False):
            if(track[ind]>0):
                track[ind]-=1
                done=True
            else:
                for i in range(ind,len(d)):
                    track[i]=d[i][1]
                ind-=1
        factor=1
        for i in range(0,len(d)):
            factor*=d[i][0]**track[i]
        f.append(factor)
    return f
def factors(n,d=0):
    """
    Returns an array of the factors of n, in ascending order.
    Do not use for things like sigma(n),
    which are significantly faster to call directly.
    """
    return sorted(unsorted_factors(n,d))
#d(n): number of diviors
def d(n,d=0):
    if(d==0):
        d=div(n)
    prod=1
    for f in d:
        prod*=f[1]+1
    return prod
def height(n,d=0):
    "Sum of prime exponents in the factorization of n"
    if(d==0):
        d=div(n)
    out=0
    for i in f:
        out+=i[1]
    return out
def sig(k,n,d=0):
    "Sigma_k(n): sum of the kth powers of the divisors of n"
    if(d==0):
        d=div(n)
    if(n==0):
        return 0
    prod=1
    for i in d:
        s=0
        for j in range (0,i[1]+1):
            s+=i[0]**(j*k)
        prod*=s
    return prod
def phi(n, d=0,fact=False):
    """Number of relatively prime integers <=n
    If fact=true, returns [phi(n),div(phi(n))]"""
    if(d==0):
        d=div(n)
    prod=1
    factorization=[]
    for f in d:
        prod*=f[0]-1
        prod*=f[0]**(f[1]-1)
        if fact:
            factorization=divmult(factorization,div(f[0]-1))
            factorization=divmult(factorization,[[f[0],f[1]-1]])
    if fact:
        return [prod,factorization]
    return prod
#squarefree part of n
def squarefree(n,d=0):
    if(d==0):
        d=div(n)
    prod=1
    for f in d:
        if(f[1]%2==1):
            prod*=f[0]
    return prod
#Mobius function: look it up on wikipedia or something
def mobius(n,d=0):
    if(d==0):
        d=div(n)
    if(squarefree(n,d)==n):
        return 1-2*(len(d)%2)
    else:
        return 0
#Mertens function: just the sum of Mobius.
#Give value m at input k if known to speed computation.
def mertens(n,k=0,m=0):
    if(k==0):
        sum=0
    else:
        sum=m
    for i in range(k+1,n+1):
        sum+=mobius(i)
    return sum
#Fibonacci with F_0=0
fibonacci=[0,1]
def fib(n):
    "Linear time, not log(n)"
    for i in range(len(fibonacci),n+1):
        fibonacci.append(fibonacci[-1]+fibonacci[-2])
    return fibonacci[n]
#Element-wise product of two lists up to k elements if desired
def lprod(l,m,k=-1):
    if(k==-1):
        k=min(len(l),len(m))
    tot=0
    for i in range(0,k):
        tot+=l[i]*m[i]
    return tot
#Generalized fibonacci: given initial terms and a recursion rule
#(e.g., [2,1]=2f(n-2)+f(n-1)), what is the nth term? s is the index
#of the first term in init.
def gfib(init,rec,n,s=0):
    l=init[:]
    for i in range(len(init),n+1-s):
        l.append(lprod(l[i-len(init):i],rec))
    return l[n]
def abundance(n,d=0):
    """Amount by which the
    the sum of factors
    of n exceeds 2n"""
    if(d==0):
        d=div(n)
    return sig(1,n,d)-2*n
def aliquot(n,d=0):
    """Sum of the proper divisors of n"""
    if(d==0):
        d=div(n)
    return sig(1,n,d)-n
def pi(n):
    "Number of primes <= n"
    return len(sieve(n+1))
def factorial(n):
    "Works on half-integers too"
    if(n==0):
        return 1
    elif(n==0.5):
        return math.sqrt(math.pi)/2
    elif(n<1):
        return 0
    else:
        return n*factorial(n-1)
def choose(n,k):
    "nCk"
    k=min(k,n-k)
    if k<0:
        return 0
    m=k
    p=1
    d=2
    while m:
        p*=n-m+1
        while d<=k and p%d==0:
            p//=d
            d+=1
        m-=1
    return p
def matmult(A,B):
    """Matrix multiplication of A and B.
    Returns -1 if dimensions don't check out"""
    #
    #
    #
    # UNFINISHED
    #
    #
    #
    #
    return -1
def digits(n,b=10):
    "Array of the digits of n in base b"
    if(n==0):
        return [0]
    dig=[]
    while(n>0):
        dig.insert(0,n%b)
        n=n//b
    return dig
def antidigits(d,b=10):
    """Convert a list of
    digits to an integer"""
    num=0
    for i in range(0,len(d)):
        num*=b
        num+=d[i]
    return num
#Number of 0s in base b 
def zeros(n,b=10):
    '''Number of 0s
    in base b'''
    if(n==0):
        return 1
    zs=0
    while(n%b==0):
        n//=b
        zs+=1
    return zs
#fill in code here - should give the smallest integer not in a list
def absent(n):
    return 0
def isum(f,s=0,e=10000):
    tot=0
    for i in range(s,e):
        tot+=f(i)
    return tot
def iprod(f,s=0,e=10000):
    tot=1
    for i in range(s,e):
        tot*=f(i)
    return tot
def continuedfraction(efrac):
    "Input a list and get [numerator,denominator] of the result"
    val=[efrac[-1],1]
    efrac=efrac[:-1]
    while(len(efrac)>0):
        val=[val[0]*efrac[-1]+val[1],val[0]]
        efrac=efrac[:-1]
    return val
def cFracRep(x,l=10):
    ind=0
    out=[]
    while(ind<l):
        out.append(int(x))
        x=x%1
        if(abs(x)>0.0000001):
            x=1/x
            ind+=1
        else:
            ind=l
    if(out[-1]==1) and len(out)>1:
        out=out[:-2]+[out[-2]+1]
    return out
def rationalCFrac(a,b):
    out=[]
    while(b!=0):
        out.append(a//b)
        a=a%b
        a=a+b
        b=a-b
        a=a-b
    return out
def dif(s):
    "First differences of a list s"
    out=[]
    for i in range(0,len(s)-1):
        out.append(s[i+1]-s[i])
    return out
def smallFactor(n):
    """Smallest prime factor of n"""
    if(n<2):
        return -1
    for p in primes:
        if n%p==0:
            return p
    t=primes[-1]+2
    while(t*t<=n):
        if n%t==0:
            return 
        t+=2
    return n
def prod(*a):
    "Product of the inputs"
    if(len(a)==1):
        a=a[0]
    p=1
    for e in a:
        p*=e
    return p
def modpow(a,e,n):
    "a^e mod n, computed efficiently"
    if(e==0):
        return 1
    init=a
    for d in bin(e)[3:]:
        a*=a
        if d=='1':
            a*=init
        a=a%n
    return a

