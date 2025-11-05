# extract code from tutorials snippets
# codeExtract file.md file.py

import re
regstart = re.compile(bytes('```python', encoding='utf-8'))
regend = re.compile(bytes('```', encoding='utf-8'))

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        print("codeExtract file1.md file2.py")
        exit(1)
     
    file1 = sys.argv[1]
    file2 = sys.argv[2]
    print("Reading %s..."%file1)
    f1 = open("%s"%file1, 'rb')
    lines = f1.readlines()
    f1.close()

    print("Writing %s..."%file2)
    f2 = open("%s"%file2, 'wb')
    
    write = False
    for l in lines:
        if regstart.search(l) is not None: 
            write = True
            f2.write(bytes("\n", encoding='utf-8'))
        elif regend.search(l) is not None: write = False
        elif write:
            f2.write(l)
    f2.close()
