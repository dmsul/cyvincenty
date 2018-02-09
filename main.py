from geopy.distance import vincenty as v1
from cyvincenty import vincenty as v2


if __name__ == '__main__':
    p1 = (34, -118)
    p2 = (35, -117)
    new = v2(p1[0], p1[1], p2[0], p2[1])
    old = v1(p1, p2).km
    print(new / old)
