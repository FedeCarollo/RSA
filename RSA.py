import random

class RSA:
    __generatingPrimes = (None, None)

    __private = None
    __public = None
    __mod = None

    def __init__(self) -> None:
        pass

    def GetPublic(self):
        return self.__public

    def SavePublicKey(self, encode, mod) -> None:
        self.__public, self.__mod = encode, mod


    def GenKeyPair(self, seed: int=3893489):
        p, q = self.GenRandomPrimes(seed, 1024)
        N = p*q

        e = self.FindEncodingExponent(p, q)
        d = self.GetInverse(e, (p-1)*(q-1))


        self.__generatingPrimes = (p, q)
        self.__private = d
        self.__public = (e, N)
        self.__mod = N

    def GetInverse(self, x, N):
        (d, _, b) = RSA.ExtendedEuclid(N, x)

        if d != 1:  #non accade
            return None
        
        return b % N
    

    def EncodeMessage(self, message):
        return RSA.FastModExp(message, self.__public[0], self.__public[1])
    
    def DecodeMessage(self, encoded):
        return RSA.FastModExp(encoded, self.__private, self.__mod)

   
    @staticmethod
    def GCD(x, y) -> int:
        if y == 0:
            return x
        
        return RSA.GCD(y, x%y)
        

    #seleziona e tc: e ≡ 1 mod 
    @staticmethod
    def FindEncodingExponent(p, q):
        mod = (p-1)*(q-1) 

        for e in range(3, mod, 2):
            if(RSA.GCD(e, mod) == 1):
                return e
            
        return None #non può accadere per costruzione


    @staticmethod
    def GenRandomPrime(seed, n=128) -> int:
        while(True):
            n_str = ""

            for _ in range(n-1):
                n_str = n_str + ("1" if random.random() >= 0.5 else "0")

            n_str+= "1" #genera sempre numeri dispari
            r = int(n_str, 2)

            if RSA.IsPrime(r):
                return r
        return None
    
    @staticmethod
    def GenRandomPrimes(seed, n=128):
        random.seed(seed)

        while(True):
            a, b = RSA.GenRandomPrime(seed), RSA.GenRandomPrime(seed)

            if a != b:
                return (a, b)
            
    @staticmethod
    def IsPrime(r: int, k:int = 20) -> bool:
        z = [random.randint(1, r) for _ in range(k)]

        for t in z:                 #per ogni sample di z
            if(RSA.GCD(t, r) > 1):      #testimone banale
                return False
            """     caso non funzionante per i numeri di Carmichael
            if(fastModExp(z, r-1, r) != 1):     #viola il piccolo teorema di fermat
                return False
            """

            divisors = []
            d = 1
            while((r-1)//d % 2 == 0):
                divisors.append((r-1)//d)
                d*=2

            for div in divisors:
                if(RSA.FastModExp(t, div, r) != 1):     #se esiste una radice non banale di 1 mod r allora r non è primo
                    return False
            
        return True     #probabilmente è primo (1/2)^k possibilità di aver sbagliato
    
    @staticmethod
    def ExtendedEuclid(x, y) -> tuple[int, int, int]:
        if(y == 0):
            return (x, 1, 0)

        (d, a, b) = RSA.ExtendedEuclid(y, x % y)
        return (d, b, a - (x // y) * b)

    @staticmethod
    def FastModExp(x, exp, N):      #result is x^exp mod N
        if exp==0:
            return 1

        a = RSA.FastModExp(x, exp//2, N)

        if(exp%2==0):
            return a**2 % N
        else:
            return x*a**2 % N