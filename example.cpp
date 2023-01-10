#include "editdis.hpp"

int main(int argc, const char** argv) {
    int dis = editdis::edit_distance("浙B82D60", "闽B82060");
    std::cout << dis << std::endl;  // 2

    // support wstring
    dis = editdis::edit_distance(L"浙B82D60", L"闽B82060");
    std::cout << dis << std::endl;  // 2

    return 0;
}