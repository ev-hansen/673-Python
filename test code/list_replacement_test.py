test_list = [[[1, 2, 3],
              [4, 5, 6]],
             [[7, 8, 9],
              [10, 11, 12]]]

new_list = [[[[k if k > 5 else float("NaN")] for k in test_list[i][j]] 
             for j in range(len(test_list[i]))] 
            for i in range(len(test_list))]

for i in range(len(new_list)):
    for j in range(len(new_list[i])):
        print(f"j: {new_list[j]}")
    print()
