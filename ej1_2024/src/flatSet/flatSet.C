#include <messageStream.H>
#include <List.H>

template <typename T>
class flatSet : public Foam::List<T>
{
private:
    Foam::List<T> elements;

public:
    flatSet() : Foam::List<T>() {};
    flatSet(Foam::List<T> l)
    {
        forAll(l, i)
        {
            append(l[i]);
        }
    }

    bool find(T e)
    {
        return elements.found(e);
    }

    void append(T e)
    {
        if (!find(e))
        {
            elements.append(e);
        }
    }

    flatSet unionSet(flatSet<T> s)
    {
        flatSet<T> result = *this;
        forAll(s.elements, i)
        {
            result.append(s.elements[i]);
        }
        return result;
    }

    flatSet interSet(flatSet<T> s)
    {
        flatSet<T> result;
        forAll(s.elements, i)
        {
            if (find(s.elements[i]))
            {
                result.append(s.elements[i]);
            }
        }
        return result;
    }
    friend Foam::Ostream& operator<<(Foam::Ostream& os, const flatSet<T>& fs)
    {
        os << fs.elements;
        return os;
    }
};

int main()
{
    Foam::List<int> l(3, 0);
    Foam::List<int> l2(3, 35);
    l.append(1);
    l.append(42);
    flatSet<int> s(l);
    s.append(35);
    s.append(42);
    s.append(36);
    flatSet<int> s2(l2);
    flatSet<int> s3 = s.interSet(s2);
    Foam::Info << s << Foam::endl;
    Foam::Info << s2 << Foam::endl;
    Foam::Info << s.unionSet(s2) << Foam::endl;
    Foam::Info << s3 << Foam::endl;
}
