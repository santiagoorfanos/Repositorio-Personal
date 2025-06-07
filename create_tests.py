r_i = "1"
r_e = "5"
iso = "500"
ninst = "5"
for i in range(1,21):
    for j in range(1,21):
        m = 5*i
        n = 5*j
        values = "1500 " * n + "100 " * n
        value_lines = [values]*int(ninst)
        with open(f"size_tests/size_test_{m}_{n}.txt","w") as file:
            lines = [r_i,r_e,str(m),str(n),iso,ninst] + value_lines
            
            file.writelines([line + '\n' for line in lines])
