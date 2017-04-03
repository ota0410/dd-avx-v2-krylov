#coding:utf-8

print("γの値を入力してください\n")
gamma = gets.chomp.to_f

#if !(gamma.is_a?(Integer)||gamma.is_a?(Float))
#  print("int型ではありません終了します\n")
#  exit(1)
#end

print("つぎに行列のサイズを入力してください\n")
n = gets.chomp.to_i

#if !(n.is_a?(Integer))
#  print("int型ではありません終了します\n")
#  exit(1)
#end

nnz = 3*n-2
file = File.open("gamma.mtx","w")

file.print("%%MatrixMarket matrix coordinate real general\n")
file.print("#{n}\s#{n}\s#{nnz}\n")
file.print("1 1 2\n")
file.print("2 1 1\n")

count = 1
while count < n-1 
  file.print("#{count}\s#{count+1}\s#{gamma}\n")
  file.print("#{count+1}\s#{count+1}\s2\n")
  file.print("#{count+2}\s#{count+1}\s1\n")
  count += 1
end

file.print("#{n-1}\s#{n}\s#{gamma}\n")
file.print("#{n}\s#{n}\s2\n")

file.close

file2 = File.open("t_gamma.mtx","w")
file2.print("%%MatrixMarket matrix coordinate real general\n")
file2.print("#{n}\s#{n}\s#{nnz}\n")
file2.print("1 1 2\n")
file2.print("2 1 #{gamma}\n")

count = 1
while count < n-1
  file2.print("#{count}\s#{count+1}\s1\n")
  file2.print("#{count+1}\s#{count+1}\s2\n")
  file2.print("#{count+2}\s#{count+1}\s#{gamma}\n")
  count += 1
end
file2.print("#{n-1}\s#{n}\s1\n")
file2.print("#{n}\s#{n}\s2\n")

file2.close
