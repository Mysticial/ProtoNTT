/* TestFramework.cpp
 * 
 * Author           : Alexander J. Yee
 * Date Created     : 11/06/2014
 * Last Modified    : 11/15/2014
 * 
 */

#include <iostream>
#include <string>
#include <chrono>
#include "TestFramework.h"
namespace ProtoNTT{
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  These functions are copy/pasted from y-cruncher.
std::string tostr_u_commas(uint64_t x){
    //  Prints out x with comma separators.

    std::string str = std::to_string(x);
    std::string out;

    const char* ptr = str.c_str();
    size_t len = str.size();

    size_t commas = (len + 2) / 3 - 1;
    size_t shift = len - commas * 3;

    while (1){
        char ch = *ptr++;
        if (ch == '\0')
            break;
        if (shift == 0){
            out += ',';
            shift = 3;
        }
        out += ch;
        shift--;
    }

    return out;
}
std::string tostr_u_bytes(uint64_t bytes){
    //  Prints out bytes in one of the following forms:
    //  0.xx suffix
    //  x.xx suffix
    //  xx.x suffix
    //   xxx suffix

    const char* BYTE_NAMES[] = {
        " bytes",
        " KiB",
        " MiB",
        " GiB",
        " TiB",
        " PiB",
        " EiB",
    };

    std::string out;
    if (bytes < 1000){
        out += std::to_string(bytes);
        out += BYTE_NAMES[0];
        return out;
    }

    int suffix_index = 1;
    while (bytes >= 1024000){
        bytes >>= 10;
        suffix_index++;
    }

    bytes *= 1000;
    bytes >>= 10;

    //  .xxx or (1.00)
    if (bytes < 1000){
        bytes += 5;
        bytes /= 10;

        if (bytes == 100){
            out += "1.00";
            out += BYTE_NAMES[suffix_index];
            return out;
        }

        out += "0.";
        out += std::to_string(bytes);
        out += BYTE_NAMES[suffix_index];
        return out;
    }

    //  x.xx or (10.0)
    if (bytes < 10000){
        bytes += 5;
        bytes /= 10;

        if (bytes == 1000){
            out += "10.0";
            out += BYTE_NAMES[suffix_index];
            return out;
        }

        out += std::to_string(bytes / 100);
        bytes %= 100;
        out += ".";
        if (bytes >= 10){
            out += std::to_string(bytes);
        }else{
            out += "0";
            out += std::to_string(bytes);
        }
        out += BYTE_NAMES[suffix_index];
        return out;
    }

    //  xx.x or (0.98)
    if (bytes < 100000){
        bytes += 50;
        bytes /= 100;

        if (bytes == 1000){
            out += " 100";
            out += BYTE_NAMES[suffix_index];
            return out;
        }

        out += std::to_string(bytes / 10);
        bytes %= 10;
        out += ".";
        out += std::to_string(bytes);
        out += BYTE_NAMES[suffix_index];
        return out;
    }

    //  xxx or (1.00)
    {
        bytes += 500;
        bytes /= 1000;

        if (bytes == 1000){
            out += "0.98";
            out += BYTE_NAMES[suffix_index + 1];
            return out;
        }

        out += " ";
        out += std::to_string(bytes);
        out += BYTE_NAMES[suffix_index];
        return out;
    }
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  Printing
void pause(){
    std::cout << "Press ENTER to continue . . .";
    std::string out;
    std::getline(std::cin,out);
}
void print_commas(uint64_t x){
    std::cout << tostr_u_commas(x);
}
uint64_t bits_to_decimal(uint64_t bits){
    const double RATIO = 0.30102999566398119521373889472449302676818988146211;
    return (uint64_t)(bits * RATIO);
}
void print_bits(uint64_t bits){
    std::cout << tostr_u_commas(bits) << " bits";
    std::cout << " ( " << tostr_u_commas(bits_to_decimal(bits)) << " decimal digits )" << std::endl;
}
void print_words(uint64_t words){
    std::cout << tostr_u_commas(words) << " words";
    std::cout << " ( " << tostr_u_commas(bits_to_decimal(words * 64)) << " decimal digits )" << std::endl;
}
void print_bytes(uint64_t bytes){
    std::cout << tostr_u_commas(bytes) << " bytes";
    std::cout << " ( " << tostr_u_bytes(bytes) << " )" << std::endl;
}
double wall_clock(){
    auto ratio_object = std::chrono::high_resolution_clock::period();
    double ratio = (double)ratio_object.num / ratio_object.den;
    return std::chrono::high_resolution_clock::now().time_since_epoch().count() * ratio;
}
void print(const uint64_t* A,size_t L){
    if (L == 0){
        std::cout << "{}" << std::endl;
        return;
    }

    std::cout << "{";

    L--;
    for (size_t c = 0; c < L; c++){
        std::cout << (A[c]) << ",";
    }
    std::cout << A[L] << "}" << std::endl;
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
}
