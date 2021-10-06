#ifndef ERROR_FUNCTION_H
#define ERROR_FUNCTION_H

#include "Connection.h"

/**
 * Base error function class to weigth the face non-orthogonality and its 
 * and contribution to the global function to optimize
 *
 * @param values Container whose values are summed.
 * @return sum of `values`, or 0.0 if `values` is empty.
 */
class Error_Function
{
    public:
    
        Error_Function() {};
        
        /**
         * Calculate the contribution of this face to the global
         * error function
         *
         * @param the face orthogonality
         * @return the contribution of the face
         */
        virtual double get_value(double a) {return 0.;}
        
        /**
         * Calculate the contribution of this face to the global
         * error function derivative
         *
         * @param the face orthogonality
         * @return the derivative of the face contribution
         */
        virtual double get_derivative(double a) {return 0.;}
};


// Orthogonality error function
/**
 * Define the error as 1-ortho
 *
 * @param No parameter
 * @return 1-orthogonality
 */
class Id_Function: public Error_Function
{
    public:
        Id_Function() {};
        double get_value(double ortho) {
            return 1-ortho;
        }
        double get_derivative(double ortho) {
            return -1.;
        }
};

/**
 * Penalize high non-orthogonality using a power function
 *
 * @param the penalization power n
 * @return (1-orthogonality)^n
 */
class Power_Function: public Error_Function
{
    public:
        double penalizing_power;
        Power_Function(double power) {
            if (power < 1) {
                std::invalid_argument("Penalizing power in power weighting face \
                error funcion cannot be null or negative");
            }
            penalizing_power = power;
        }
        
        double get_value(double ortho) { 
            if (ortho <= 0.) {
                return std::log(-1); //NaN
            }
            else {
                return std::pow(1-ortho, penalizing_power);
            }
        }
        double get_derivative(double ortho) {
            return - penalizing_power * std::pow(1-ortho, penalizing_power-1);
        }
};

class Inverse_Function: public Error_Function
{
    public:
        double penalizing_power;
        Inverse_Function(double power) {
            if (power <= 0.) {
                std::invalid_argument("Penalizing power in inverse weighting \
                face error funcion cannot be null or negative");
            }
            penalizing_power = power;
        }
        
        double get_value(double ortho) {
            if (ortho <= 0.) {
                return std::log(-1); //NaN
            }
            else {
                return 1/std::pow(ortho, penalizing_power) - 1;
            }
        }
        double get_derivative(double ortho) {
            return - penalizing_power * std::pow(ortho, penalizing_power+1);
        }
};

class Log_Function: public Error_Function
{
    public:
        Log_Function() {};
        double get_value(double ortho) {
            return -std::log(ortho);
        }
        double get_derivative(double ortho) {
            return -1/ortho;
        }
};

class Tan_Function: public Error_Function
{
    public:
        Tan_Function() {};
        double get_value(double ortho) {
            return -std::tan(ortho);
        }
        double get_derivative(double ortho) {
            double cos = std::cos(ortho);
            return -1/(cos*cos);
        }
};

#endif // ERROR_FUNCTION_H
