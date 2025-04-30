#include <Field.H>

class myVectorField : public Foam::Field<Foam::vector>
{
public:
    using Foam::Field<Foam::vector>::Field;

    
    Foam::Field<Foam::scalar> operator&(const myVectorField& vecField)
    {
       Foam::Field<Foam::scalar> component_0(this->size());
       Foam::Field<Foam::scalar> component_1(this->size());
       Foam::Field<Foam::scalar> component_2(this->size());
       component_0 = this->component(0) * vecField.component(0);
       component_1 = this->component(1) * vecField.component(1);
       component_2 = this->component(2) * vecField.component(2);

       return component_0 + component_1 + component_2;
    }
    
    Foam::Field<Foam::scalar> mag()
    {
        Foam::Field<Foam::scalar> prodPunto = *this&*this;
        return sqrt(prodPunto);
    }

    myVectorField operator*(const myVectorField& vecField)
    {
        myVectorField crossProduct(this->size());
        Foam::Field<Foam::scalar> comp0 = this->component(1) * vecField.component(2) - 
                                          this->component(2) * vecField.component(1);
        Foam::Field<Foam::scalar> comp1 = this->component(2) * vecField.component(0) - 
                                          this->component(0) * vecField.component(2);
        Foam::Field<Foam::scalar> comp2 = this->component(0) * vecField.component(1) - 
                                          this->component(1) * vecField.component(0);

        crossProduct.replace(0, comp0);
        crossProduct.replace(1, comp1);
        crossProduct.replace(2, comp2);

        return crossProduct;
    }
};

int main()
{
    Foam::vector v1(1,2,3);
    Foam::vector v2(4,5,6);
    Foam::vector v3(7,8,9);
    Foam::vector v4(10,11,12);
    Foam::List<Foam::vector> vecList(4);
    vecList[0] = v1;
    vecList[1] = v2;
    vecList[2] = v3;
    vecList[3] = v2;

    Foam::List<Foam::vector> vecList2(4);
    vecList2[0] = v3;
    vecList2[1] = v4;
    vecList2[2] = v1;
    vecList2[3] = v2;

    myVectorField vecField(vecList);
    myVectorField vecField2(vecList2);

    Foam::Field<Foam::scalar> prodPunto = vecField&vecField2;
    myVectorField prodCruz = vecField*vecField2;
    Foam::Field<Foam::scalar>  prodPunto2 = prodCruz&vecField2;

    Foam::Info << "Vector Field 1: " << vecField << Foam::endl;
    Foam::Info << "Vector Field 2: " << vecField2 << Foam::endl;
    Foam::Info << "Suma: " << vecField+vecField2 << Foam::endl;
    Foam::Info << "Resta: " << vecField-vecField2 << Foam::endl;
    Foam::Info << "10xField1: " << 10*vecField << Foam::endl;
    Foam::Info << "Producto punto: " << prodPunto << Foam::endl;
    Foam::Info << "Magnitud VectorField 1: " << vecField.mag() << Foam::endl;
    Foam::Info << "Producto cruz: " << prodCruz << Foam::endl;
    Foam::Info << "Magnitud prodCruz: " << prodCruz.mag() << Foam::endl;
    Foam::Info << "Producto punto 2: " << prodPunto2 << Foam::endl;
}
