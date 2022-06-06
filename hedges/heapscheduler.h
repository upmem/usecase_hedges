/* usage:
HeapScheduler<Doub,CargoClass> myheap;
myheap.push(time, cargo);
timeval = myheap.pop(cargoval);  // cargo is returned in cargoval
// or
timeval = myheap.pop();
cargoval = myheap.lastcargo;
// or, if no cargo, can just do
HeapScheduler<> myheap;  // (equivalent to <Doub,void*>, actually)
myheap.push(time);
timeval = myheap.pop();
*/
typedef int32_t heap_score_type;
typedef uint32_t heap_ptr_type;

#define FLOAT_TO_FP(i) (heap_score_type)((float)(i) * ((float)(1 << DECODER_QUANT_FRAC_BITS)))
#define FP_TO_FLOAT(i) (float)((float)(i) / ((float)(1 << DECODER_QUANT_FRAC_BITS)))
#define HEAP_SCORE_MAX_VAL_FLOAT 32767


template <class T = Doub, class U = void *>
struct HeapScheduler
{
public:
	static const Int defaultps = 1100000; // initial heap size
	T bigval;
	U lastcargo;
	Int ps, ks;
	NRvector<T> ar; // times
	NRvector<U> br; // "cargo"

	HeapScheduler() : bigval(HEAP_SCORE_MAX_VAL_FLOAT), ps(defaultps), ks(0), ar(ps, bigval), br(ps) {}
	void push(T time, U cargo = U(NULL))
	{
		// lengthen list, add to end, sift up
		// pushes a time and cargo onto the heap
		int32_t ark, armo;
		Int k, mo;
		if (ks == ps)
			resizear(2 * ps);

		k = ks++;
		ar[k] = time;
		br[k] = cargo;
		mo = (k - 1) / 2;
		ark = FLOAT_TO_FP(ar[k]);
		armo = FLOAT_TO_FP(ar[mo]);
		while (k > 0 && armo > ark)
		{
			SWAP(ar[k], ar[mo]); // swap with mother
			SWAP(br[k], br[mo]);
			k = mo;
			mo = (k - 1) / 2;
			ark = FLOAT_TO_FP(ar[k]);
			armo = FLOAT_TO_FP(ar[mo]);
		}
	}
	uint64_t size()
	{
		return (uint64_t)(ks);
	}

	T pop() { return pop(lastcargo); } // if no argument, return cargo in HeapScheduler::lastcargo
	T pop(U &cargo)
	{ // return top of heap, move last to top, shorten list, sift down
		// pops the next (in order) time and its cargo from the heap
		// returns numeric_limits::max() time, and U() cargo, when heap is empty
		Int k = 0, rdau, ldau, mindau;
		T ans = ar[0];
		U cans = br[0];

		int32_t ark, arr, arl, mindauu;
		if ((ks--) > 0)
		{
			ar[0] = ar[ks];
			br[0] = br[ks];
			ar[ks] = bigval;
			br[ks] = U();
			while ((ldau = 2 * k + 1) < ks)
			{
				rdau = ldau + 1; // might be a bigval, but that is OK
				ark = FLOAT_TO_FP(ar[k]);
				arr = FLOAT_TO_FP(ar[rdau]);
				arl = FLOAT_TO_FP(ar[ldau]);
				mindau = (arl < arr ? ldau : rdau);
				mindauu = FLOAT_TO_FP( ar[mindau]);
				if (ark > mindauu)
				{
					// printf("[SWAP] %lu %lu \n", k, mindau);
					SWAP(ar[k], ar[mindau]); // swap with smaller of two daughters
					SWAP(br[k], br[mindau]);
				}
				else
					break;
				k = mindau;
			}
		}
		cargo = cans;
		return ans;
	}
	void resizear(Int newps)
	{							// only used internally
		ar.resize(newps, true); // resize preserving contents
		br.resize(newps, true);
		for (int i = ks; i < newps; i++)
			ar[i] = bigval;
		ps = newps;
	}
	void rewind()
	{ // zero out the heap w/o changing its size in memory
		ks = 0;
		ar.assign(ps, bigval);
		br.resize(ps);
		// for (int i = 0; i < ps; i++) {
		//	ar[i] = bigval;
		//	br[i] = U(NULL);
		// }
	}
	void reinit()
	{ // zero out the heap and give back memory
		ks = 0;
		ps = defaultps;
		ar.assign(ps, bigval);
		br.resize(ps);
	}

	void printheap()
	{ // only used for debugging
		// for (int i = 0; i < ks; i++)
		for (int i = 0; i < 100; i++)
			printf("HEAP [%d]  %.5f  %u\n", i, ar[i], (unsigned)(br[i]));
	}
};
