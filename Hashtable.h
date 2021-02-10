class Hashtable {
    Hashtable *one;
    Hashtable *zero;
    int first;
 public:
    Hashtable(unsigned long);
    ~Hashtable();
    void free_table(Hashtable*);
    void store_individual(unsigned long *, unsigned long);
    unsigned long check_individual(unsigned long *, unsigned long);
    unsigned long delete_individual(unsigned long *, unsigned long);
    void reset(unsigned long);
    unsigned long nrElem;
    unsigned long nrIndividuals;

    Hashtable** listCurrent;
    unsigned long* indexCurrent;
    long iterator(unsigned long *, unsigned long);
    long iteratorInit(unsigned long*, unsigned long);
    void iteratorFree(unsigned long);

    Hashtable** secondListCurrent;
    unsigned long* secondIndexCurrent;
    long iteratorSecond(unsigned long *, unsigned long);
    long iteratorSecondInit(unsigned long*, unsigned long);
};

