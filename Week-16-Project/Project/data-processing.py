f = open('data-L-vary.txt')

dm = 1
L = 30

data = f.readlines()
for i in range(len(data)):
     data[i] = data[i].split()
     if (data[i] == []):
          print('{0:>10d}{1:>20.10f}{2:>20.10f}{3:>20.10f}{4:>20.10f}'.format(L,T_max,m_max,chi_max,C_max))
          dm = 1
          L += 2
     else:
          for j in range(4):
              data[i][j] = float(data[i][j])
          if (abs(data[i][1] - 0.5) < dm):
              dm = abs(data[i][1] - 0.5)
              T_max = data[i][0]
              m_max = data[i][1]
              chi_max = data[i][2]
              C_max = data[i][3]

f.close()
print()
