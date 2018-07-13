include ../Make.inc

OBJS = send_buffer.o recv_buffer.o

default: libBuffer.a


libBuffer.a: $(OBJS)
	ar rcs $@ $^

cleanall:
	@rm -f *.o *~ *.a *.mod
