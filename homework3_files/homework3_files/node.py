#https://www.tutorialspoint.com/python/python_nodes.htm
#https://www.geeksforgeeks.org/self-in-python-class/

class Node:
    def __init__(self, id):
        self.id = id
        self.children = None

    def add_child(self, node, distance):
        if not self.children:
            self.children = {}
        self.children[node] = distance
