class doubleList{
	double* object;
	unsigned long numbers;
	doubleList* next;
	unsigned long index;
	public:
		doubleList();
		doubleList(double);
		doubleList(const doubleList&);
		void setObject(double);
		int setObject(double,unsigned long);
		double getObject();
		double getObject(unsigned long);
		unsigned long getNumbers(unsigned long);
		int addElement(double);
		void print();
		//ostream& operator<<(doubleList*);
		~doubleList();
		double operator[](int);
		unsigned long operator()(int);
		int operator==(doubleList&);
		int operator!=(doubleList&);
		unsigned long getIndex(double);
		int contains(double); 
		unsigned long getSize();
};
