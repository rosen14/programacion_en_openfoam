#include <iostream>

template<typename T>
class lista {
protected:
    struct Nodo {
        T dato;
        Nodo* siguiente;
    };

    Nodo* cabeza;

public:
    lista() : cabeza(nullptr){}
    ~lista() { //destructor
        Nodo* temp;
        while (cabeza != nullptr) {
            temp = cabeza;
            cabeza = cabeza->siguiente;
            delete temp; // elimino desde la cabeza
        }
    }
    void append(T dato) { //inserto al final de la lista
        Nodo* nuevo = new Nodo{dato, nullptr};
        if (cabeza == nullptr) { // si la lista está vacia, el nuevo nodo se convierte en la cabeza
            cabeza = nuevo;
        } else {
            Nodo* temp = cabeza;
            while (temp->siguiente != nullptr) {
                temp = temp->siguiente;
            }
            temp->siguiente = nuevo;
        }
    }
    // Dado que esta función no es miembro de la clase uso friend para que tenga acceso a los miembros
    // publicos y protegidos
    friend std::ostream& operator<<(std::ostream& os, const lista<T>& l) {
        Nodo* temp = l.cabeza;
        while (temp != nullptr) {
            os << temp->dato << " ";
            temp = temp->siguiente;
        }
        return os;
    }
};  

template <typename T>
class pila : public lista<T> { // va a ser una pila LIFO
public:
    void push(T dato) {
        this->append(dato);
    }
    T top() { // Devuelve el último elemento ingresado sin eliminarlo
        if (this->cabeza == nullptr) { // si no hay elementos en la pila
            throw std::out_of_range("Pila vacia");
        }
        typename lista<T>::Nodo* temp = this->cabeza;
        while (temp->siguiente != nullptr) { // Avanza hasta el último nodo
            temp = temp->siguiente;
        }
        return temp->dato; // Devuelve el dato del último nodo
    }
    bool empty() { // Devuelve true si la pila está vacía
        return this->cabeza == nullptr;
    }
    T pop() { // Elimina el último elemento ingresado y devuelve el valor
        if (this->cabeza == nullptr) { // si no hay elementos en la pila
            throw std::out_of_range("Pila vacia");
        }
    
        typename lista<T>::Nodo* temp = this->cabeza;
        T dato;
    
        if (this->cabeza->siguiente == nullptr) { // Si solo hay un nodo
            dato = this->cabeza->dato;
            delete this->cabeza;
            this->cabeza = nullptr;
        } else { // Si hay más de un nodo
            while (temp->siguiente->siguiente != nullptr) { // Avanza hasta el penúltimo nodo
                temp = temp->siguiente;
            }
            dato = temp->siguiente->dato; // Guarda el dato del último nodo
            delete temp->siguiente;      // Elimina el último nodo
            temp->siguiente = nullptr;  // Actualiza el penúltimo nodo para que sea el último
        }
    
        return dato;
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
};

int main() {
    pila<int> p;
    p.push(3);
    p.push(5);
    p.push(1);
    p.push(4);
    p.push(-4);
    p.push(0);
    std::cout << "Pila antes de ordenar: " << p << std::endl;

    p.sort();

    std::cout << "Pila después de ordenar: " << p << std::endl;

    return 0;
}
