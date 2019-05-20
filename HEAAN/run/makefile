# This sample should be run after make libHEAAN.a file in lib folder

# All Target
all: clean TestHEAAN

# Tool invocations
TestHEAAN:
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C++ Linker'
	g++ test.cpp -o TestHEAAN -L/usr/local/lib ../lib/libHEAAN.a -lntl -lgmp -lm -std=c++11 -O2 -I/../src -pthread
	#g++ toy.cpp -o toy -L/usr/local/lib ../lib/libHEAAN.a -lntl -lgmp -lm -std=c++11 -O2 -I/../src -pthread
	g++ ToyRingMultiplier.cpp toy.cpp -o toy -L/usr/local/lib -lntl -lgmp -lm -std=c++11 -O2 -I/. -pthread
	@echo 'Finished building target: $@'
	@echo ' '

# Other Targets
clean:
	rm -rf TestHEAAN
	rm -rf toy
