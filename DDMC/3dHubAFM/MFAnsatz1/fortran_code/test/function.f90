program test_function
implicit none
real*8 x,y,temp
real*8 add,subt
x=1
y=2
print *, add(x,y)
print *, subt(x,y)
end program

function subt(x,y)
implicit none
real*8 subt,add,x,y
subt=add(x,y)-x
return
end function



function add(x,y)
implicit none
real*8 add,x,y
add=x+y
return
end function

