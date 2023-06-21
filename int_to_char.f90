subroutine int_to_char(number,str3)
implicit none
integer*4, intent(in)::number
integer*4 :: num
integer::aux, count, i
character(1),dimension(10)::str = (/'1','2','3','4','5','6','7','8','9','0'/)
character(9), intent(out)::str3
	num = number
	str3 = ''

	if (num >= 100000000) then
		count = 9
	else if (num >= 10000000) then
		count = 8
	else if (num >= 1000000) then
		count = 7
	else if (num >= 100000) then
		count = 6
	else if (num >= 10000) then
		count = 5
	else if (num >= 1000) then
		count = 4
	else if (num >= 100) then
		count = 3
	else if (num >= 10) then
		count = 2
	else
		count = 1
	end if
	do i=1,6
		if (i <= count) then
			aux = mod(num,10)
			if (aux == 0) then
				str3 = str(10)//str3
			else
				str3 = str(aux)//str3
			end if
			num = num/10
		else
			str3 = str(10)//str3
		end if
	end do
return
end subroutine int_to_char
