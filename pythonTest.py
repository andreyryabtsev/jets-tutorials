def test(a):
    print(a + 1)

def test2(b):
    b.increment()
    print(b)





class TestObject:
    def __init__(self, val):
        self.value = val
    def __str__(self):
        return str(self.value)
    def increment(self):
        self.value = self.value + 1
    def __add__(self, val):
        self.value = self.value + val
        return self



a = TestObject(3)
test(a)
test2(a)




