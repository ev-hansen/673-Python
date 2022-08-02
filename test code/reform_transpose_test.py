import numpy


def reform(given_list):
    arr = numpy.array(given_list)
    reformed = arr.squeeze()
    reformed_list = reformed.tolist()
    return reformed_list


def transpose(given_list):
    arr = numpy.array(given_list)
    transposed = arr.T
    transposed_list = transposed.tolist()
    return transposed_list


test_list_1 = [[1], [2, 3], 4]
ref_test_1 = reform(test_list_1)
trans_test_1 = transpose(test_list_1)
print(f"list: {test_list_1}\n reformed: {ref_test_1}\n"
      f"transposed: {trans_test_1}")

test_list_2 = [[1, 2], [3, 4], [5, 6]]
ref_test_2 = reform(test_list_2)
trans_test_2 = transpose(test_list_2)
print(f"list: {test_list_2}\n reformed: {ref_test_2}\n"
      f"transposed: {trans_test_2}")

test_list_3 = [[1, 2], [3], [4]]
ref_test_3 = reform(test_list_3)
trans_test_3 = transpose(test_list_3)
print(f"list: {test_list_3}\n reformed: {ref_test_3}\n"
      f"transposed: {trans_test_3}")

test_list_4 = [[1], [2], [3]]
ref_test_4 = reform(test_list_4)
trans_test_4 = transpose(test_list_4)
print(f"list: {test_list_4}\n reformed: {ref_test_4}\n"
      f"transposed: {trans_test_4}")
