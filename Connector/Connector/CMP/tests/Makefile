include ../Make.inc

OBJ1 = test_shared_ptr.o
OBJ2 = test_vector_view.o
OBJ3 = test_pending_message_container.o

default:	test_pending_message_container.exe  test_shared_ptr.exe  test_vector_view.exe	

test_pending_message_container.exe: $(OBJ3)
	$(CXX) -o $@ $< $(LINK) 

test_shared_ptr.exe: $(OBJ1)
	$(CC) -o $@ $<   $(LINK)

test_vector_view.exe: $(OBJ2)
	$(F77) -o $@ $<  $(LINK)

cleanall:
	@rm -f *.o *~ *.exe
