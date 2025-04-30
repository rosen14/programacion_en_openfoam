#include <List.H>

template<typename T>
class pila : private Foam::List<T> // No quiero exponer las funcionalidades de List a los usuarios de pila
{
public:
    using Foam::List<T>::List; // Heredar constructores de Foam::List<T>

    void push(T item)
    {
        this->append(item); // Usar append de Foam::List<T>
    }

    T top()
    {
        if (this->empty())
        {
            Foam::FatalError << "Pila vacía" << Foam::endl;
        }
        return this->last(); // método de UList<T>
    }
    T pop()
    {
        T item = this->top();
        this->pop_back(); // método de List<T>
        return item;
    }

    void sort() {
        pila<T> tempStack; // Pila temporal para almacenar elementos ordenados
    
        // Mientras la pila original no esté vacía
        while (!this->empty()) {
            T current = this->pop(); // Saca el elemento superior de la pila original
            
            // Mueve elementos de tempStack de vuelta a la pila original si son mayores que el actual
            while (!tempStack.empty() && tempStack.top() > current) {
                this->push(tempStack.pop());
                
            }
    
            // Inserta el elemento actual en la pila temporal
            tempStack.push(current);
        }
    
        // Mueve los elementos de la pila temporal de vuelta a la pila original
        while (!tempStack.empty()) {
            this->push(tempStack.pop());
        }
    }
    friend Foam::Ostream& operator<<(Foam::Ostream& os, const pila<T>& p)
    {
        for (int i = 0; i < p.size(); ++i)
        {
            os << p[i] << " ";
        }
        return os;
    }
  };

int main()
{
    pila<int> p;
    p.push(3);
    p.push(2);
    p.push(1);
    Foam::Info << p << Foam::endl;
    Foam::Info << "Pop: " << p.pop() << Foam::endl;
    Foam::Info << p << Foam::endl;
    Foam::Info << "Pop: " << p.pop() << Foam::endl;
    Foam::Info << p << Foam::endl;
    p.push(5);
    p.push(4);
    p.push(-9);
    Foam::Info << p << Foam::endl;
    p.sort();
    Foam::Info << "Sorted: " << p << Foam::endl;
}
