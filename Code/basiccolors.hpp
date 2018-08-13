#include <vector>
using namespace std;
vector<int> brewer_sequential = {5, 44, 67, 6, 48, 72, 8, 51, 76, 10, 55, 81, 11, 58, 85, 13, 61, 90, 16, 65, 94, 18, 68, 98, 20, 71, 103, 22, 75, 107, 24, 78, 111, 26, 81, 115, 29, 84, 119, 31, 87, 123, 33, 90, 127, 35, 94, 131, 38, 97, 135, 40, 100, 139, 43, 103, 143, 45, 106, 146, 47, 109, 150, 50, 112, 154, 52, 114, 157, 55, 117, 161, 58, 120, 164, 60, 123, 167, 63, 126, 171, 65, 128, 174, 68, 131, 177, 71, 134, 180, 73, 136, 183, 76, 139, 186, 79, 141, 189, 81, 144, 191, 84, 147, 194, 87, 149, 197, 90, 151, 199, 93, 154, 202, 95, 156, 204, 98, 158, 206, 101, 161, 209, 104, 163, 211, 107, 165, 213, 110, 167, 215, 113, 170, 217, 116, 172, 219, 119, 174, 220, 122, 176, 222, 125, 178, 224, 128, 180, 226, 130, 182, 227, 133, 184, 229, 136, 186, 230, 139, 188, 232, 142, 190, 233, 145, 192, 234, 148, 194, 235, 150, 195, 237, 153, 197, 238, 156, 199, 239, 159, 201, 240, 161, 202, 241, 164, 204, 242, 167, 206, 243, 169, 207, 244, 172, 209, 245, 175, 211, 245, 177, 212, 246, 180, 214, 247, 182, 215, 248, 185, 217, 248, 187, 218, 249, 190, 220, 249, 192, 221, 250, 195, 222, 251, 197, 224, 251, 199, 225, 252, 202, 227, 252, 204, 228, 252, 206, 229, 253, 209, 231, 253, 211, 232, 253, 213, 233, 254, 215, 234, 254, 218, 236, 254, 220, 237, 254, 222, 238, 255, 224, 239, 255, 226, 240, 255, 228, 242, 255, 230, 243, 255, 232, 244, 255, 234, 245, 255, 236, 246, 255, 238, 247, 255, 240, 248, 255, 242, 249, 255, 244, 250, 255, 246, 251, 255, 248, 252, 255};
vector<int> brewer_diverging = {5, 44, 67, 8, 51, 76, 11, 58, 85, 16, 65, 94, 20, 71, 103, 24, 78, 111, 29, 84, 119, 33, 90, 127, 38, 97, 135, 43, 103, 143, 47, 109, 150, 52, 114, 157, 58, 120, 164, 63, 126, 171, 68, 131, 177, 73, 136, 183, 79, 141, 189, 84, 147, 194, 90, 151, 199, 95, 156, 204, 101, 161, 209, 107, 165, 213, 113, 170, 217, 119, 174, 220, 125, 178, 224, 130, 182, 227, 136, 186, 230, 142, 190, 233, 148, 194, 235, 153, 197, 238, 159, 201, 240, 164, 204, 242, 169, 207, 244, 175, 211, 245, 180, 214, 247, 185, 217, 248, 190, 220, 249, 195, 222, 251, 199, 225, 252, 204, 228, 252, 209, 231, 253, 213, 233, 254, 218, 236, 254, 222, 238, 255, 226, 240, 255, 230, 243, 255, 234, 245, 255, 238, 247, 255, 242, 249, 255, 246, 251, 255, 238, 255, 232, 231, 255, 223, 224, 254, 215, 218, 253, 207, 212, 252, 199, 207, 251, 192, 201, 249, 185, 196, 247, 179, 191, 246, 172, 186, 243, 166, 181, 241, 160, 176, 239, 154, 171, 236, 148, 167, 234, 142, 162, 231, 136, 158, 228, 131, 153, 225, 125, 149, 222, 120, 144, 218, 115, 140, 215, 110, 136, 211, 105, 132, 207, 100, 127, 204, 95, 123, 200, 90, 119, 195, 86, 115, 191, 81, 111, 187, 77, 108, 182, 73, 104, 177, 69, 100, 173, 65, 96, 168, 61, 92, 162, 57, 88, 157, 53, 84, 152, 49, 80, 146, 46, 76, 141, 42, 72, 135, 38, 68, 129, 35, 64, 123, 31, 60, 117, 28, 56, 110, 25, 52, 104, 21, 48, 97, 18, 44, 91, 15, 40, 84, 12, 35, 77, 9, 31, 70, 6, 27, 63, 5, 23, 55, 3, 19, 48, 2};
vector<int> brewer_qualitative = {237, 92, 130, 238, 94, 125, 240, 97, 121, 241, 100, 116, 242, 103, 111, 242, 105, 106, 243, 108, 100, 244, 111, 93, 244, 114, 87, 244, 117, 79, 245, 120, 71, 245, 123, 61, 245, 127, 49, 244, 130, 33, 242, 134, 29, 238, 138, 34, 235, 143, 38, 232, 146, 42, 230, 150, 45, 227, 153, 48, 225, 157, 51, 223, 160, 53, 221, 163, 56, 219, 166, 58, 217, 168, 60, 216, 171, 62, 214, 174, 64, 212, 176, 66, 211, 179, 69, 209, 182, 71, 207, 184, 73, 206, 187, 75, 204, 189, 77, 202, 192, 80, 200, 194, 82, 198, 196, 84, 195, 197, 85, 190, 197, 85, 185, 196, 84, 180, 196, 84, 175, 195, 83, 169, 195, 83, 164, 194, 83, 158, 194, 82, 153, 193, 82, 147, 193, 81, 140, 192, 81, 133, 192, 80, 126, 191, 80, 118, 191, 79, 109, 191, 78, 99, 190, 78, 87, 190, 77, 77, 189, 81, 78, 187, 94, 79, 185, 103, 80, 183, 110, 81, 181, 116, 81, 179, 120, 82, 177, 124, 83, 175, 127, 83, 173, 130, 84, 171, 132, 84, 169, 134, 84, 167, 136, 85, 165, 137, 85, 163, 138, 85, 161, 139, 85, 160, 140, 86, 158, 140, 86, 156, 141, 86, 154, 141, 86, 153, 141, 86, 151, 141, 86, 149, 141, 86, 147, 141, 86, 146, 141, 86, 144, 141, 86, 142, 141, 86, 141, 140, 85, 139, 140, 84, 138, 140, 83, 136, 140, 82, 135, 140, 81, 133, 139, 81, 131, 139, 80, 130, 139, 79, 128, 139, 78, 127, 139, 77, 125, 139, 76, 124, 138, 75, 122, 138, 74, 120, 139, 74, 119, 139, 73, 117, 139, 72, 116, 139, 71, 114, 140, 70, 112, 141, 69, 111, 142, 67, 109, 143};
vector<int> plsequential_lightness = {0, 0, 0, 12, 1, 0, 21, 3, 1, 27, 5, 1, 32, 6, 1, 37, 8, 2, 40, 9, 3, 44, 11, 3, 47, 13, 4, 50, 14, 5, 53, 16, 6, 56, 17, 7, 59, 19, 8, 62, 21, 9, 65, 22, 11, 69, 24, 12, 72, 26, 14, 75, 28, 15, 78, 29, 17, 81, 31, 19, 84, 33, 20, 87, 35, 22, 90, 37, 24, 93, 39, 26, 97, 41, 28, 100, 43, 30, 103, 45, 31, 106, 47, 33, 109, 49, 35, 112, 51, 37, 115, 53, 39, 118, 55, 42, 121, 57, 44, 124, 60, 46, 127, 62, 48, 130, 64, 50, 133, 66, 53, 136, 69, 55, 139, 71, 57, 142, 73, 59, 145, 76, 62, 148, 78, 64, 151, 80, 67, 154, 83, 69, 157, 85, 72, 160, 88, 74, 163, 90, 77, 165, 93, 79, 168, 95, 82, 171, 98, 85, 174, 100, 87, 176, 103, 90, 179, 106, 93, 182, 108, 96, 184, 111, 98, 187, 114, 101, 190, 116, 104, 192, 119, 107, 195, 122, 110, 197, 125, 113, 200, 128, 116, 202, 130, 119, 204, 133, 122, 207, 136, 125, 209, 139, 128, 211, 142, 131, 214, 145, 134, 216, 148, 138, 218, 151, 141, 220, 154, 144, 222, 157, 147, 224, 160, 151, 226, 163, 154, 228, 166, 157, 230, 169, 161, 232, 173, 164, 233, 176, 168, 235, 179, 171, 237, 182, 175, 238, 185, 178, 240, 189, 182, 241, 192, 185, 243, 195, 189, 244, 199, 193, 245, 202, 196, 247, 205, 200, 248, 209, 204, 249, 212, 208, 250, 216, 211, 251, 219, 215, 252, 223, 219, 252, 226, 223, 253, 230, 227, 254, 233, 231, 254, 237, 235, 254, 240, 239, 255, 244, 243, 255, 248, 247, 255, 251, 251, 255, 255, 255};
vector<int> plsequential_saturation = {206, 76, 29, 206, 77, 31, 205, 77, 33, 204, 78, 35, 204, 78, 37, 203, 79, 39, 203, 79, 40, 202, 80, 42, 201, 80, 43, 201, 81, 45, 200, 81, 46, 200, 82, 48, 199, 82, 49, 198, 83, 50, 198, 83, 52, 197, 84, 53, 196, 84, 54, 196, 85, 56, 195, 85, 57, 195, 85, 58, 194, 86, 59, 193, 86, 60, 193, 87, 61, 192, 87, 62, 191, 88, 63, 191, 88, 64, 190, 89, 65, 189, 89, 66, 189, 90, 67, 188, 90, 68, 187, 91, 69, 186, 91, 70, 186, 91, 71, 185, 92, 72, 184, 92, 73, 184, 93, 74, 183, 93, 75, 182, 94, 76, 181, 94, 77, 181, 95, 77, 180, 95, 78, 179, 95, 79, 178, 96, 80, 178, 96, 81, 177, 97, 82, 176, 97, 82, 175, 98, 83, 174, 98, 84, 174, 98, 85, 173, 99, 86, 172, 99, 86, 171, 100, 87, 170, 100, 88, 169, 101, 89, 169, 101, 89, 168, 101, 90, 167, 102, 91, 166, 102, 92, 165, 103, 92, 164, 103, 93, 163, 103, 94, 162, 104, 95, 161, 104, 95, 161, 105, 96, 160, 105, 97, 159, 106, 97, 158, 106, 98, 157, 106, 99, 156, 107, 99, 155, 107, 100, 154, 108, 101, 153, 108, 101, 152, 108, 102, 151, 109, 103, 150, 109, 103, 149, 110, 104, 148, 110, 105, 147, 110, 105, 145, 111, 106, 144, 111, 107, 143, 112, 107, 142, 112, 108, 141, 112, 109, 140, 113, 109, 139, 113, 110, 138, 114, 110, 136, 114, 111, 135, 114, 112, 134, 115, 112, 133, 115, 113, 131, 115, 114, 130, 116, 114, 129, 116, 115, 127, 117, 115, 126, 117, 116, 125, 117, 117, 123, 118, 117, 122, 118, 118, 120, 119, 118, 119, 119, 119};
vector<int> plsequential_rainbow = {0, 0, 0, 9, 2, 4, 16, 4, 8, 22, 6, 13, 26, 8, 18, 30, 11, 22, 33, 13, 26, 35, 15, 30, 37, 17, 34, 38, 19, 37, 40, 21, 41, 41, 23, 44, 42, 26, 48, 42, 28, 51, 42, 31, 54, 42, 34, 57, 42, 37, 60, 41, 40, 63, 40, 43, 66, 38, 46, 68, 36, 49, 71, 34, 52, 73, 31, 55, 75, 29, 59, 76, 25, 62, 78, 22, 65, 79, 18, 68, 80, 15, 71, 81, 12, 74, 81, 9, 77, 82, 9, 80, 82, 10, 83, 82, 13, 85, 82, 18, 88, 81, 23, 91, 81, 29, 93, 80, 35, 95, 80, 41, 98, 79, 47, 100, 78, 53, 102, 78, 60, 105, 77, 66, 107, 77, 72, 109, 76, 79, 111, 76, 85, 112, 76, 91, 114, 76, 98, 116, 77, 104, 118, 78, 110, 119, 79, 116, 121, 81, 122, 123, 84, 127, 124, 86, 133, 126, 89, 139, 127, 93, 144, 129, 97, 149, 131, 101, 154, 132, 105, 159, 134, 110, 164, 136, 115, 168, 137, 120, 172, 139, 126, 176, 141, 131, 180, 143, 137, 184, 145, 142, 187, 147, 148, 190, 150, 154, 192, 152, 159, 195, 154, 165, 197, 157, 170, 199, 160, 175, 200, 163, 180, 201, 166, 185, 203, 169, 190, 203, 172, 194, 204, 176, 199, 205, 179, 203, 205, 183, 206, 205, 186, 210, 206, 190, 213, 206, 193, 216, 206, 197, 219, 207, 200, 222, 207, 204, 224, 208, 207, 226, 209, 211, 228, 210, 214, 230, 212, 218, 232, 213, 221, 233, 215, 224, 235, 218, 227, 236, 220, 230, 238, 223, 233, 239, 226, 236, 241, 230, 239, 242, 234, 242, 244, 238, 245, 246, 242, 247, 248, 246, 250, 250, 251, 253, 252, 255, 255, 255};
vector<int> plsequential_blackbody = {0, 0, 0, 13, 1, 2, 23, 2, 4, 29, 3, 6, 34, 4, 8, 39, 5, 9, 43, 6, 10, 46, 8, 11, 49, 10, 12, 52, 11, 13, 55, 13, 13, 58, 15, 14, 61, 17, 14, 64, 19, 14, 67, 21, 14, 70, 23, 14, 72, 25, 14, 75, 28, 13, 78, 30, 13, 81, 32, 12, 83, 34, 12, 86, 37, 11, 89, 39, 10, 91, 41, 9, 94, 44, 8, 96, 46, 7, 99, 49, 6, 102, 51, 5, 104, 54, 4, 107, 56, 4, 109, 59, 3, 112, 61, 3, 114, 64, 3, 117, 66, 3, 119, 69, 4, 122, 71, 5, 124, 74, 6, 127, 76, 7, 129, 79, 9, 132, 82, 12, 134, 84, 14, 137, 87, 17, 139, 89, 20, 142, 92, 23, 144, 95, 25, 147, 97, 28, 149, 100, 31, 152, 103, 35, 154, 105, 38, 156, 108, 41, 159, 111, 44, 161, 113, 47, 164, 116, 51, 166, 119, 54, 169, 122, 58, 171, 124, 61, 173, 127, 65, 176, 130, 68, 178, 132, 72, 181, 135, 75, 183, 138, 79, 185, 141, 83, 188, 144, 87, 190, 146, 91, 192, 149, 95, 195, 152, 99, 197, 155, 103, 199, 157, 107, 202, 160, 111, 204, 163, 115, 206, 166, 119, 209, 169, 123, 211, 172, 128, 213, 174, 132, 216, 177, 136, 218, 180, 141, 220, 183, 145, 222, 186, 150, 225, 189, 154, 227, 192, 159, 229, 194, 164, 231, 197, 169, 233, 200, 173, 236, 203, 178, 238, 206, 183, 240, 209, 188, 242, 212, 193, 244, 215, 199, 246, 218, 204, 248, 221, 209, 250, 224, 215, 251, 227, 220, 253, 230, 226, 254, 233, 231, 255, 236, 237, 255, 240, 242, 255, 243, 247, 255, 247, 251, 255, 251, 253, 255, 255, 255};
vector<int> pldiverging_lightness = {37, 7, 1, 44, 10, 2, 50, 13, 3, 56, 16, 5, 62, 19, 7, 68, 22, 9, 74, 26, 12, 80, 29, 15, 86, 33, 18, 92, 36, 22, 98, 40, 25, 104, 44, 29, 110, 47, 33, 116, 51, 36, 122, 55, 40, 127, 59, 45, 133, 64, 49, 139, 68, 53, 145, 72, 58, 150, 77, 62, 156, 82, 67, 161, 86, 72, 167, 91, 77, 172, 96, 82, 177, 101, 87, 182, 106, 92, 188, 111, 97, 192, 116, 103, 197, 121, 109, 202, 127, 114, 207, 132, 120, 211, 138, 126, 215, 143, 132, 219, 149, 138, 223, 155, 144, 227, 161, 151, 231, 167, 157, 234, 173, 164, 237, 179, 171, 240, 185, 177, 243, 191, 184, 246, 198, 191, 248, 204, 198, 250, 211, 206, 252, 217, 213, 253, 224, 220, 254, 231, 228, 255, 237, 235, 255, 244, 243, 255, 251, 251, 252, 252, 255, 245, 246, 254, 239, 241, 254, 232, 235, 253, 226, 229, 252, 219, 224, 250, 213, 218, 249, 207, 212, 247, 200, 207, 245, 194, 201, 243, 188, 196, 241, 182, 190, 238, 176, 185, 235, 170, 180, 232, 164, 174, 229, 158, 169, 226, 152, 164, 223, 146, 158, 219, 141, 153, 215, 135, 148, 211, 130, 143, 207, 124, 138, 203, 119, 133, 199, 113, 128, 195, 108, 123, 190, 103, 118, 185, 98, 113, 181, 93, 108, 176, 87, 103, 171, 83, 99, 166, 78, 94, 160, 73, 89, 155, 68, 85, 150, 64, 80, 144, 59, 76, 139, 55, 71, 133, 50, 67, 127, 46, 62, 121, 42, 58, 115, 38, 54, 109, 34, 50, 103, 30, 46, 97, 26, 42, 91, 23, 38, 85, 19, 34, 79, 16, 30, 72, 13, 26, 66, 10, 23, 60, 7, 19, 53, 5, 15, 45};
vector<int> pldiverging_saturation = {206, 76, 29, 205, 77, 33, 204, 78, 37, 203, 79, 40, 201, 80, 43, 200, 81, 46, 199, 82, 49, 198, 83, 52, 196, 84, 54, 195, 85, 57, 194, 86, 59, 193, 87, 61, 191, 88, 63, 190, 89, 65, 189, 90, 67, 187, 91, 69, 186, 91, 71, 184, 92, 73, 183, 93, 75, 181, 94, 77, 180, 95, 78, 178, 96, 80, 177, 97, 82, 175, 98, 83, 174, 98, 85, 172, 99, 86, 170, 100, 88, 169, 101, 89, 167, 102, 91, 165, 103, 92, 163, 103, 94, 161, 104, 95, 160, 105, 97, 158, 106, 98, 156, 107, 99, 154, 108, 101, 152, 108, 102, 150, 109, 103, 148, 110, 105, 145, 111, 106, 143, 112, 107, 141, 112, 109, 139, 113, 110, 136, 114, 111, 134, 115, 112, 131, 115, 114, 129, 116, 115, 126, 117, 116, 123, 118, 117, 120, 119, 118, 119, 119, 120, 118, 119, 123, 118, 119, 126, 117, 119, 128, 116, 118, 131, 116, 118, 133, 115, 118, 136, 115, 118, 138, 114, 118, 141, 113, 118, 143, 113, 118, 146, 112, 117, 148, 111, 117, 151, 111, 117, 153, 110, 117, 155, 109, 117, 158, 108, 117, 160, 108, 117, 163, 107, 116, 165, 106, 116, 167, 105, 116, 170, 104, 116, 172, 103, 116, 174, 102, 116, 177, 101, 115, 179, 100, 115, 182, 99, 115, 184, 98, 115, 186, 97, 115, 189, 96, 115, 191, 95, 114, 193, 94, 114, 196, 93, 114, 198, 91, 114, 201, 90, 114, 203, 89, 113, 206, 87, 113, 208, 86, 113, 210, 84, 113, 213, 83, 112, 215, 81, 112, 218, 79, 112, 220, 77, 112, 223, 75, 111, 225, 73, 111, 228, 71, 111, 231, 69, 111, 233, 66, 110, 236, 64, 110, 238, 61, 110, 241};
vector<int> plqualitative_hue = {195, 104, 124, 194, 105, 119, 194, 106, 114, 193, 107, 109, 191, 108, 104, 190, 110, 99, 188, 111, 94, 187, 112, 89, 185, 114, 84, 183, 115, 78, 181, 117, 73, 178, 118, 67, 176, 120, 62, 173, 121, 56, 171, 123, 50, 168, 124, 44, 165, 126, 38, 162, 127, 32, 158, 128, 26, 155, 130, 19, 151, 131, 12, 148, 133, 6, 144, 134, 2, 140, 135, 1, 136, 136, 2, 131, 138, 5, 127, 139, 11, 122, 140, 18, 117, 141, 24, 112, 142, 31, 106, 143, 37, 100, 144, 43, 94, 145, 49, 87, 146, 54, 80, 147, 60, 72, 148, 66, 63, 149, 71, 53, 149, 76, 41, 150, 82, 23, 151, 87, 0, 151, 92, 0, 152, 97, 0, 152, 102, 0, 153, 107, 0, 153, 112, 0, 153, 116, 0, 153, 121, 0, 154, 125, 0, 154, 130, 0, 154, 134, 0, 154, 139, 0, 153, 143, 0, 153, 147, 0, 153, 151, 0, 152, 155, 0, 152, 159, 0, 151, 162, 0, 150, 166, 0, 149, 169, 0, 149, 172, 0, 148, 175, 0, 146, 178, 0, 145, 181, 0, 144, 184, 0, 142, 186, 29, 141, 188, 49, 139, 190, 63, 138, 192, 75, 136, 193, 85, 134, 194, 94, 132, 195, 103, 130, 196, 111, 128, 197, 119, 126, 197, 126, 124, 197, 132, 122, 197, 138, 120, 196, 144, 118, 196, 150, 116, 195, 155, 114, 194, 160, 112, 192, 164, 111, 190, 168, 109, 189, 172, 107, 186, 176, 106, 184, 179, 105, 182, 182, 103, 179, 184, 102, 176, 187, 101, 173, 189, 101, 169, 190, 100, 166, 192, 100, 162, 193, 100, 159, 194, 100, 155, 195, 100, 151, 196, 100, 146, 196, 101, 142, 196, 101, 138, 196, 102, 133, 196, 103, 129};
vector<int> cubehelix = {0, 0, 0, 5, 1, 4, 9, 2, 9, 13, 4, 14, 17, 5, 19, 20, 7, 25, 22, 9, 30, 24, 12, 36, 26, 14, 42, 27, 17, 47, 27, 20, 53, 27, 23, 58, 27, 27, 63, 27, 31, 68, 26, 35, 72, 24, 39, 75, 23, 43, 78, 22, 48, 81, 20, 53, 83, 19, 58, 84, 17, 63, 85, 16, 68, 85, 15, 73, 85, 14, 78, 84, 14, 82, 83, 14, 87, 81, 14, 92, 78, 15, 96, 75, 17, 100, 72, 19, 104, 69, 22, 108, 65, 25, 111, 62, 29, 114, 58, 33, 117, 54, 39, 119, 51, 44, 121, 48, 51, 123, 45, 57, 125, 42, 65, 126, 40, 72, 126, 38, 81, 127, 37, 89, 127, 36, 98, 127, 37, 107, 127, 37, 116, 126, 39, 125, 125, 41, 134, 125, 44, 142, 124, 48, 151, 123, 52, 160, 122, 58, 168, 121, 63, 175, 120, 70, 183, 120, 77, 189, 119, 85, 196, 119, 93, 201, 119, 101, 206, 119, 110, 211, 119, 119, 214, 120, 129, 217, 121, 138, 220, 123, 148, 222, 124, 157, 223, 126, 166, 223, 129, 175, 223, 132, 184, 223, 135, 193, 222, 138, 201, 220, 142, 208, 218, 146, 215, 216, 150, 222, 214, 154, 227, 211, 159, 233, 209, 164, 237, 206, 169, 241, 203, 174, 244, 201, 179, 247, 199, 184, 249, 197, 189, 250, 195, 195, 251, 194, 199, 252, 193, 204, 252, 192, 209, 251, 192, 214, 250, 193, 218, 249, 194, 222, 248, 195, 226, 247, 197, 229, 245, 200, 233, 244, 203, 236, 243, 206, 239, 242, 210, 241, 241, 214, 243, 241, 219, 245, 241, 224, 247, 241, 229, 249, 242, 234, 250, 244, 240, 252, 246, 245, 253, 248, 250, 254, 251, 255, 255, 255};
vector<int> moreland = {180, 4, 38, 184, 20, 41, 188, 31, 44, 192, 39, 46, 195, 47, 49, 199, 54, 52, 202, 60, 55, 206, 67, 59, 209, 73, 62, 212, 78, 65, 215, 84, 68, 218, 89, 72, 221, 95, 75, 223, 100, 79, 226, 105, 82, 228, 110, 86, 231, 115, 89, 233, 120, 93, 235, 124, 97, 237, 129, 100, 238, 133, 104, 240, 138, 108, 241, 142, 112, 243, 147, 116, 244, 151, 120, 245, 155, 124, 245, 159, 128, 246, 163, 132, 247, 167, 136, 247, 170, 140, 247, 174, 144, 247, 177, 148, 247, 181, 152, 247, 184, 156, 247, 187, 160, 246, 190, 164, 246, 193, 168, 245, 196, 172, 244, 199, 176, 243, 201, 180, 241, 204, 184, 240, 206, 188, 238, 208, 192, 236, 210, 196, 234, 212, 200, 232, 214, 204, 230, 216, 208, 227, 217, 212, 225, 219, 215, 222, 220, 219, 219, 221, 222, 217, 220, 225, 214, 220, 228, 211, 219, 231, 209, 218, 233, 206, 217, 236, 203, 216, 238, 200, 215, 240, 197, 214, 243, 194, 212, 244, 190, 211, 246, 187, 209, 248, 184, 207, 249, 180, 205, 250, 177, 203, 252, 174, 201, 253, 170, 199, 253, 167, 196, 254, 163, 194, 254, 160, 191, 255, 156, 188, 255, 153, 186, 255, 149, 183, 255, 146, 180, 254, 142, 176, 254, 138, 173, 253, 135, 170, 252, 131, 166, 251, 128, 163, 250, 124, 159, 249, 121, 156, 247, 117, 152, 246, 114, 148, 244, 110, 144, 242, 107, 140, 240, 103, 136, 238, 100, 132, 235, 97, 128, 233, 93, 124, 230, 90, 120, 227, 87, 116, 224, 84, 111, 221, 80, 107, 218, 77, 103, 215, 74, 98, 211, 71, 94, 208, 68, 89, 204, 65, 85, 200, 62, 81, 196, 59, 76, 192};
vector<int> mcnames = {0, 0, 0, 0, 8, 0, 0, 16, 2, 0, 23, 5, 0, 29, 10, 0, 34, 16, 0, 38, 23, 0, 40, 31, 0, 41, 39, 0, 41, 49, 0, 40, 58, 0, 37, 68, 0, 34, 78, 0, 29, 87, 0, 23, 96, 0, 17, 104, 3, 10, 111, 11, 3, 117, 21, 0, 121, 33, 0, 124, 45, 0, 125, 58, 0, 125, 72, 0, 123, 86, 0, 120, 100, 0, 115, 114, 0, 109, 128, 0, 101, 141, 0, 93, 153, 0, 83, 164, 0, 73, 174, 0, 62, 183, 5, 52, 189, 17, 41, 194, 30, 30, 197, 45, 21, 199, 60, 12, 198, 77, 4, 195, 93, 0, 191, 111, 0, 185, 128, 0, 177, 145, 0, 168, 161, 0, 157, 177, 0, 146, 192, 0, 134, 205, 1, 122, 217, 9, 109, 227, 19, 96, 236, 31, 84, 242, 44, 72, 247, 59, 62, 250, 75, 52, 250, 92, 44, 249, 109, 37, 245, 127, 32, 240, 145, 29, 233, 162, 28, 225, 180, 29, 216, 196, 31, 205, 212, 36, 194, 226, 43, 182, 239, 51, 170, 250, 62, 158, 255, 73, 146, 255, 86, 135, 255, 100, 125, 255, 115, 115, 255, 131, 107, 255, 147, 101, 255, 163, 96, 255, 179, 92, 255, 194, 91, 255, 209, 91, 255, 223, 93, 249, 235, 96, 240, 247, 101, 231, 255, 108, 222, 255, 116, 213, 255, 126, 205, 255, 136, 197, 255, 147, 190, 255, 159, 184, 255, 171, 179, 255, 183, 175, 255, 195, 173, 255, 206, 172, 255, 217, 173, 255, 227, 175, 255, 236, 178, 255, 243, 183, 255, 250, 188, 253, 255, 195, 249, 255, 202, 247, 255, 210, 245, 255, 218, 244, 255, 227, 245, 255, 235, 247, 255, 242, 250, 255, 249, 255, 255, 255};
