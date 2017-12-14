import os


def readFile(filePath):
    try:
        file_object = open(filePath, "r", encoding="utf-8")
        all_the_text = []
        line = str.lower(str.strip(file_object.readline()))
        while line:
            all_the_text.append(line)
            line = str.lower(str.strip(file_object.readline()))
    except IOError:
        print("Error: 没有找到文件或读取文件失败")
    finally:
        if file_object is not None:
            file_object.close()
    return all_the_text


def readAllFile(filePath):
    file_object = open(filePath, "r", encoding="utf-8")
    try:
        all_the_text = file_object.read()
    finally:
        file_object.close()

    return all_the_text
