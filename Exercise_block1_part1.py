### Exercise 1 ###
def get_sphere_volume(radius):
   """Returns the volume of a sphere"""
    volume = (4/3) * 3.14 * radius ** 3
    return volume

#vol = get_sphere_volume(4)
#print (vol)

### Exercise 2 ###
def recusive_factorial(n):
    """Calculates the factorial of a number recursively"""
    if n==1:
        return n
    else:
        return n*recusive_factorial(n-1)

#fact = recusive_factorial(4)
#print (fact)

def factorial(n):
    """Calculates the factorial of a number"""
    factorial = 1
    while n > 0:
        factorial = factorial * n
        n = n-1
    return factorial

#fact = factorial(4)
#print (fact)

### Exercise 3 ###
def recursive_count_up(n,odd):
    """Counts up in the screen but when odd is true it only prints odd numbers recursively"""
    if odd == True:
        if n >= 0:
            recursive_count_up(n-1,odd)
            if n % 2 == 0:
                print (n)
    else:
        if n >= 0:
            recursive_count_up(n-1,odd)
            print (n)

#recursive_count_up(5,True)

def count_up(n,odd):
    """Counts up in the screen but when odd is true it only prints odd numbers"""
    if odd == True:
        p = 1
        while p <= n:
            print (p)
            p = p + 2
    else:
        p = 0
        while p <= n:
            print (p)
            p = p + 1

#count_up(4,True)

### Exercise 4 ###
def get_final_price(price, discount_percentage=10): # the order of the variables was incorrect
    """Return the final price after applying the discount percentage"""
    return price - ((price * discount_percentage) / 100) # the formula did not work and the variable name was incorrect

#final_price = get_final_price(40)
#print (final_price)
