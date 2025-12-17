#include <iostream>
#include <fstream>
#include <io.h>
#include <fcntl.h>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <string>

using namespace std;

// -- Helper Functions --

void printVector(vector<int> vec) {
	for (int i = 0; i < vec.size(); i++) {
		if (i != 0) {
			wcout << ", ";
		}
		wcout << vec[i];
	}
}

void printVector(vector<bool> vec) {
	for (int i = 0; i < vec.size(); i++) {
		if (i != 0) {
			wcout << ", ";
		}
		wcout << vec[i];
	}
}

template<typename T>
bool isIn2DVector(T element, vector<vector<T>> vec) {
	bool found = false;
	
	for (int i = 0; i < vec.size(); i++) {
		for (int j = 0; j < vec[i].size(); j++) {
			if (element == vec[i][j]) {
				found = true;
				break;
			}
		}
		
		if (found) {
			break;
		}
	}
	
	return found;
}

bool anyTrue(vector<bool> vec) {
	bool toReturn = false;
	for (int i = 0; i < vec.size(); i++) {
		if (vec[i]) {
			toReturn = true;
			break;
		}
	}
	
	return toReturn;
}

bool anyTrueCol(vector<vector<bool>> vec2D, int col) {
	if (vec2D.size() == 0) return false;
	if (col < 0 || col >= vec2D.size()) return false;
	
	bool toReturn = false;
	for (int i = 0; i < vec2D[0].size(); i++) {
		if (vec2D[col][i]) {
			toReturn = true;
			break;
		}
	}
	
	return toReturn;
}

vector<vector<bool>> readMaskFromFile(string filename) {
	vector<vector<bool>> mask = {};
	
	ifstream fileReader(filename);
	string line;
	
	int length;
	vector<bool> row;
	while (getline (fileReader, line)) {
		length = line.length();
		row = {};
		for (int i = 0; i < length; i++) {
			if (line[i] == '1') row.push_back(true);
			else row.push_back(false);
		}
		mask.push_back(row);
	}
	
	fileReader.close();
	return mask;
}

// -- Data Structure Classes --

template<typename T>
class Stack {
	public:
		vector<T> stack;
		int size;
		
		Stack() {
			stack = {};
			size = 0;
		}
		
		void push(T elem) {
			size++;
			
			stack.push_back(elem);
		}
		
		T peek() {
			return stack.back();
		}
		
		T pop() {
			size--;
			
			T toReturn = stack.back();
			stack.pop_back();
			return toReturn;
		}
		
		bool isInStack(T test) {
			bool found = false;
			
			for (int i = 0; i < stack.size(); i++) {
				if (test == stack[i]) {
					found = true;
					break;
				}
			}
			
			return found;
		}
};

template<typename T>
class Queue {
	public:
		vector<T> queue;
		int size;
		
		Queue() {
			queue = {};
			size = 0;
		}
		
		void enqueue(T elem) {
			size++;
			
			queue.insert(queue.begin(), elem);
		}
		
		T dequeue() {
			size--;
			
			T toReturn = queue.back();
			queue.pop_back();
			return toReturn;
		}
		
		T peek() {
			return queue.back();
		}
		
		bool isInQueue(T test) {
			bool found = false;
			
			for (int i = 0; i < queue.size(); i++) {
				if (test == queue[i]) {
					found = true;
					break;
				}
			}
			
			return found;
		}
};

class Coordinates {
	public:
		int x;
		int y;
		
		Coordinates() {
			x = 0;
			y = 0;
		}
		
		Coordinates(int num1, int num2) {
			x = num1;
			y = num2;
		}
		
		bool operator==(const Coordinates& obj) const
		{
			return x == obj.x && y == obj.y;
		}
		
		friend wostream& operator<<(wostream& wos, const Coordinates& dt);
};

wostream& operator<<(wostream& wos, const Coordinates& coor)
{
    wos << '(' << coor.x << ", " << coor.y << ')';
    return wos;
}

class TreeCellBool {
	public:
		bool element;
		vector<TreeCellBool> children;
		TreeCellBool* parent;
		
		TreeCellBool(bool elem) {
			element = elem;
			parent = NULL;
		}
		
		TreeCellBool(bool elem, TreeCellBool* p) {
			element = elem;
			parent = p;
		}
		
		void AddChild(bool elem) {
			TreeCellBool newCell(elem, this);
			children.push_back(newCell);
		}
		
		void print(int indent = 0, int lastIndent = -1) {
			if (indent == lastIndent) {
				wcout << ", ";
			} else {
				wcout << "\n";
				for (int i = 0; i < indent; i++) {
					if (i != indent - 1) {
						wcout << "  ";
					} else {
						wcout << "- ";
					}
				}
			}
			
			wcout << element;
			
			if (children.size() == 1) {
				children[0].print(indent, indent);
			} else {
				for (int i = 0; i < children.size(); i++) {
					children[i].print(indent + 1, indent);
				}
			}
		}
};

class TreeCell {
	public:
		Coordinates element;
		vector<TreeCell> children;
		TreeCell* parent;
		
		TreeCell(Coordinates elem) {
			element = elem;
			element.x = elem.x;
			element.y = elem.y;
			parent = NULL;
		}
		
		TreeCell(Coordinates elem, TreeCell* p) {
			element = elem;
			element.x = elem.x;
			element.y = elem.y;
			parent = p;
		}
		
		void AddChild(Coordinates elem) {
			TreeCell newCell(elem, this);
			children.push_back(newCell);
		}
		
		void print(int indent = 0, int lastIndent = -1) {
			if (indent == lastIndent) {
				wcout << ", ";
			} else {
				wcout << "\n";
				for (int i = 0; i < indent; i++) {
					if (i != indent - 1) {
						wcout << "  ";
					} else {
						wcout << "- ";
					}
				}
			}
			
			wcout << element;
			
			if (children.size() == 1) {
				children[0].print(indent, indent);
			} else {
				for (int i = 0; i < children.size(); i++) {
					children[i].print(indent + 1, indent);
				}
			}
		}
};

// -- Maze Class --

class Maze {
	public:
		int length;
		int width;
		vector<vector<bool>> mask;
		Coordinates start;
		Coordinates end;
		// horizontalWalls: Wether or not there is a wall above the jth cell in the ith row.
		// - dimensions: (length+1) * width
		vector<vector<bool>> horizontalWalls;
		// verticalWalls: Wether or not there is a wall to the left of the jth cell in the ith column.
		// - dimensions: (width+1) * length
		vector<vector<bool>> verticalWalls;
		// To-do: make verticalWalls indexing consistent with horizontalWalls or vice versa
		
		Maze(int s) {
			length = s;
			width = s;
			generateDefaultMask();
			generateMaze();
		}
		
		Maze(int l, int w) {
			length = l;
			width = w;
			generateDefaultMask();
			generateMaze();
		}
		
		Maze(vector<vector<bool>> m) {
			mask = m;
			if (validate_mask()) {
				length = mask.size();
				width = mask[0].size();
				generateMaze();
			} else {
				length = 5;
				width = 5;
				generateDefaultMask();
				generateMaze();
			}
		}
		
		void generateDefaultMask() {
			vector<vector<bool>> def_mask = {};
			vector<bool> temp;
			
			for (int i = 0; i < length; i++) {
				temp = {};
				for (int j = 0; j < width; i++) {
					temp.push_back(true);
				}
				def_mask.push_back(temp);
			}
			
			mask = def_mask;
			start.x = 0;
			start.y = 0;
			end.x = width - 1;
			end.y = length - 1;
		}
		
		bool validate_mask() {
			// First, ensure it has anything at all.
			int temp_length = mask.size();
			if (mask.size() == 0) {
				return false;
			}
			
			// Second, ensure that all sub-rows are the same length.
			// Also, make sure that each sub-row's length is greater than 0.
			int temp_width = mask[0].size();
			if (temp_width == 0) return false;
			for (int i = 1; i < mask.size(); i++) {
				if (mask[i].size() != temp_width) return false;
			}
			
			// Then, make sure that all true cells are connected.
			// Start by finding the first true cell.
			// Find the last true cell while we're at it.
			Coordinates temp_start;
			bool first_found = false;
			Coordinates temp_end;
			for (int i = 0; i < temp_length; i++) {
				for (int j = 0; j < temp_width; j++) {
					if (mask[i][j]) {
						temp_end.x = j;
						temp_end.y = i;
						if (!first_found) {
							temp_start.x = j;
							temp_start.y = i;
							first_found = true;
						}
					}
				}
			}
			
			Queue<Coordinates> toSearch;
			vector<vector<bool>> searched = {};
			int true_cells = 0;
			for (int i = 0; i < temp_length; i++) {
				vector<bool> temp = {};
				for (int j = 0; j < temp_width; j++) {
					temp.push_back(false);
					
					if (mask[i][j]) true_cells++;
				}
				searched.push_back(temp);
			}
			
			searched[temp_start.y][temp_start.x] = true;
			toSearch.enqueue(temp_start);
			Coordinates current;
			int x;
			int y;
			int cells_searched = 0;
			
			while (toSearch.size > 0) {
				current = toSearch.dequeue();
				x = current.x;
				y = current.y;
				
				if (x > 0 && mask[y][x - 1] && !searched[y][x - 1]) {
					Coordinates next(x - 1, y);
					toSearch.enqueue(next);
					searched[y][x - 1] = true;
				}
				
				if (y > 0 && mask[y - 1][x] && !searched[y - 1][x]) {
					Coordinates next(x, y - 1);
					toSearch.enqueue(next);
					searched[y - 1][x] = true;
				}
				
				if (x < temp_width - 1 && mask[y][x + 1] && !searched[y][x + 1]) {
					Coordinates next(x + 1, y);
					toSearch.enqueue(next);
					searched[y][x + 1] = true;
				}
				
				if (y < temp_length - 1 && mask[y + 1][x] && !searched[y + 1][x]) {
					Coordinates next(x, y + 1);
					toSearch.enqueue(next);
					searched[y + 1][x] = true;
				}
				
				cells_searched++;
			}
			
			if (cells_searched != true_cells) return false;
			
			// Finally, ensure that the mask has at least one true value on the topmost and bottommost rows...
			//        ... as well as the leftmost and rightmost columns.
			// - Just "shave" off all rows and columns at the edges that are made entirely of false values.
			while (mask.size() > 0 && !anyTrue(mask[0])) {
				mask.erase(mask.begin());
				temp_length -= 1;
				temp_start.y -= 1;
				temp_end.y -= 1;
			}
			
			while (mask.size() > 0 && !anyTrue(mask[mask.size() - 1])) {
				mask.pop_back();
				temp_length -= 1;
			}
			
			if (mask.size() == 0) return false;
			
			while (mask[0].size() > 0 && !anyTrueCol(mask, 0)) {
				for (int i = 0; i < mask.size(); i++) {
					mask[i].erase(mask[i].begin());
				}
				temp_width -= 1;
				temp_start.x -= 1;
				temp_end.x -= 1;
			}
			
			while (mask[0].size() > 0 && !anyTrueCol(mask, mask.size() - 1)) {
				for (int i = 0; i < mask.size(); i++) {
					mask[i].pop_back();
				}
				temp_width -= 1;
			}
			
			if (mask[0].size() == 0) return false;
			
			// Now that this mask is valid, assign all temp variables permanently!
			length = temp_length;
			width = temp_width;
			start = temp_start;
			end = temp_end;
			
			return true;
		}
		
		// ALGORITHM:
		// - Start with a square grid with a wall separating each adjaccent tile.
		// - Establish a start cell from where the maze is entered.
		// - Establish an end cell from where the maze is exited.
		// - Perform a depth-first search from the start cell with few twists.
		//   - Greedily search random cells out of all available candidates.
		//     - Available candidates are unsearched and are labelled as "true" in the mask.
		//   - Open the wall between cells when travelling to a new adjacent cell.
		//   - Every random range of new tiles discovered, save the current cell into a list of "split cells".
		//     - Travel back to these cells when encountering a dead end to generate a split path from the split cell.
		//     - dfs backtracking will remember these warps and "reverse" them later.
		//       - Example: If you warp from dead end (13, 2) to split cell (5, 5)...
		//                  ...you will eventually backtrack directly from (5, 5) to (13, 2).
		void generateMaze() {
			
			horizontalWalls = {};
			verticalWalls = {};
			vector<vector<bool>> cellsWithPath = {};
			Stack<Coordinates> pathTaken;
			
			int randomOffset = 0;
			
			// Start with all cells surrounded by walls.
			for (int i = 0; i < length + 1; i++) {
				horizontalWalls.push_back({});
				for (int j = 0; j < width; j++)
					horizontalWalls[i].push_back(true);
			}
			
			for (int i = 0; i < width + 1; i++) {
				verticalWalls.push_back({});
				if (i != width)
					cellsWithPath.push_back({});
				for (int j = 0; j < length; j++) {
					verticalWalls[i].push_back(true);
					if (i != width)
						cellsWithPath[i].push_back(false);
				}
			}
			
			// Remove the top wall from the start cell and the bottom wall of the end cell.
			horizontalWalls[start.y][start.x] = false;
			cellsWithPath[start.x][start.y] = true;
			horizontalWalls[end.y + 1][end.x] = false;
			
			Coordinates currentCell(start.x, start.y);
			cellsWithPath[start.x][start.y] = true;
			vector<int> validDirections;
			int pathToTake;
			
			// Keep a list of points we've travelled along to retrun to to make path splits.
			vector<Coordinates> pathSplits = {};
			
			// The distance till next path split is a random number between
			//   5 and twice the the maze size
			//   - size: midpoint between length and width
			//   strictly 5 if the maze size is 2 or 1 (not that it matters).
			int size = (length + width) / 2;
			int splitMax = size;
			int splitMin = 5;
			if (splitMin > splitMax)
				splitMax = splitMin;
			int tillNextPathSplit = (rand() % (splitMax - splitMin)) + 1 + splitMin;
			int x;
			int y;
			
			int iterations = 0;
			
			do {
				// Get a list of valid directions to travel from here.
				x = currentCell.x;
				y = currentCell.y;
				validDirections = {};
				
				// Give the up and left directions a higher chance...
				//   ...to be picked to stall going to the end.
				if (y > 0 && !cellsWithPath[x][y - 1] && mask[y - 1][x]) {
					validDirections.push_back(0);
					validDirections.push_back(0);
				}
				
				if (x > 0 && !cellsWithPath[x - 1][y] && mask[y][x - 1]) {
					validDirections.push_back(1);
					validDirections.push_back(1);
				}
			
				if (y < length - 1 && !cellsWithPath[x][y + 1] && mask[y + 1][x])
					validDirections.push_back(2);
					
				if (x < width - 1 && !cellsWithPath[x + 1][y] && mask[y][x + 1])
					validDirections.push_back(3);
				
				if (validDirections.size() > 0 && !(currentCell == end)) {
					if (tillNextPathSplit == 0) {
						// Go back to a split path and set the next counter
						pathSplits.push_back(currentCell);
						tillNextPathSplit = (rand() % (splitMax - splitMin)) + splitMin;
					} else {
						// Pick a random cell and travel there
						pathToTake = validDirections[rand() % validDirections.size()];
						// Remember the path that took us here so we can eventually backtrack.
						// This guarantees that all valid cells in the mask are a part of some path.
						pathTaken.push(currentCell);
						
						if (pathToTake == 0) {
							horizontalWalls[currentCell.y][currentCell.x] = false;
							currentCell.y = currentCell.y - 1;
						}
						
						if (pathToTake == 1) {
							verticalWalls[currentCell.x][currentCell.y] = false;
							currentCell.x = currentCell.x - 1;
						}
						
						if (pathToTake == 2) {
							horizontalWalls[currentCell.y + 1][currentCell.x] = false;
							currentCell.y = currentCell.y + 1;
						}
						
						if (pathToTake == 3) {
							verticalWalls[currentCell.x + 1][currentCell.y] = false;
							currentCell.x = currentCell.x + 1;
						}
						
						cellsWithPath[currentCell.x][currentCell.y] = true;
						tillNextPathSplit--;
					}
				} else {
					// If we found a dead end...
					// - (the end cell counts as a dead end)
					if (pathTaken.size > 0) {
						if (pathSplits.size() != 0) {
							// Generate a split path if any split path markers are left
							currentCell.x = pathSplits[pathSplits.size() - 1].x;
							currentCell.y = pathSplits[pathSplits.size() - 1].y;
							pathSplits.pop_back();
						} else {
							// Backtrack
							currentCell.x = pathTaken.peek().x;
							currentCell.y = pathTaken.peek().y;
							pathTaken.pop();
						}
					}
				}
				
				iterations++;
				// exit once we backtracked all the way back to the start.
			} while ((currentCell.x != start.x || currentCell.y != start.y));
			
			wcout << "Ended after " << iterations << " iterations." << endl;
		}
		
		// Print a list of all maze cell coordinates and their active walls.
		void exploreCells() {
			bool walls;
			for (int i = 0; i < width; i++) {
				for (int j = 0; j < length; j++) {
					walls = false;
					wcout << "Cell (" << i << "," << j << ") walls: ";
					
					if (horizontalWalls[i][j]) {
						wcout << "t";
						walls = true;
					}
					
					if (verticalWalls[j][i]) {
						if (walls) {
							wcout << ", ";
						}
						wcout << "l";
						walls = true;
					}
					
					if (verticalWalls[j + 1][i]) {
						if (walls) {
							wcout << ", ";
						}
						wcout << "r";
						walls = true;
					}
					
					if (horizontalWalls[i + 1][j]) {
						if (walls) {
							wcout << ", ";
						}
						wcout << "b";
						walls = true;
					}
					wcout << endl;
				}
			}
			wcout << endl;
		}
		
		// Base case of getPathsRecur(...)
		TreeCell getPaths() {
			TreeCell currentCell(start);
			
			vector<vector<bool>> checked = {};
			vector<bool> temp;
			for (int i = 0; i < width; i++) {
				temp = {};
				for (int j = 0; j < length; j++) {
					temp.push_back(false);
				}
				checked.push_back(temp);
			}
			
			checked[start.x][start.y] = true;
			
			getPathsRecur(currentCell, checked);
			
			return currentCell;
		}
		
		// Build and return a tree of cells representing all possible paths to take from the start cell.
		void getPathsRecur(TreeCell &currentCell, vector<vector<bool>> &checked) {
			
			Coordinates current = currentCell.element;
			Coordinates upper(current.x, current.y - 1);
			Coordinates left(current.x - 1, current.y);
			Coordinates lower(current.x, current.y + 1);
			Coordinates right(current.x + 1, current.y);
			
			if (upper.y >= 0 && !checked[upper.x][upper.y]
			&& !horizontalWalls[current.y][current.x]) {
				currentCell.AddChild(upper);
				checked[upper.x][upper.y] = true;
			}
			
			if (left.x >= 0 && !checked[left.x][left.y]
			&& !verticalWalls[current.x][current.y]) {
				currentCell.AddChild(left);
				checked[left.x][left.y] = true;
			}
			
			if (lower.y < length && !checked[lower.x][lower.y]
			&& !horizontalWalls[current.y + 1][current.x]) {
				currentCell.AddChild(lower);
				checked[lower.x][lower.y] = true;
			}
			
			if (right.x < width && !checked[right.x][right.y]
			&& !verticalWalls[current.x + 1][current.y]) {
				currentCell.AddChild(right);
				checked[right.x][right.y] = true;
			}
			
			for (int i = 0; i < currentCell.children.size(); i++) {
				getPathsRecur(currentCell.children[i], checked);
			}
		}
		
		// Find the end cell using a dfs. Print the path taken.
		vector<int> depthFirstSearch() {
			TreeCell paths = getPaths();
			TreeCell* currentCell = &paths;
			
			TreeCellBool checked(false);
			TreeCellBool* currentChecked = &checked;
			int nextChild = 0;
			
			vector<int> path = {};
			int iterations = 0;
			
			while((!(currentCell->element == end)
			&& !(nextChild == -1 && currentCell->element == paths.element))) {
				wcout << currentCell->element << ", ";
				
				while (currentChecked->children.size() < currentCell->children.size()) {
					currentChecked->AddChild(false);
				}
				
				nextChild = -1;
				for (int i = 0; i < currentChecked->children.size(); i++) {
					if (!(currentChecked->children[i].element)) {
						nextChild = i;
						break;
					}
				}
				
				if (nextChild == -1) {
					currentCell = currentCell->parent;
					currentChecked->element = true;
					currentChecked = currentChecked->parent;
					path.pop_back();
				} else {
					path.push_back(nextChild);
					currentCell = &(currentCell->children[nextChild]);
					currentChecked = &(currentChecked->children[nextChild]);
				}
				
				iterations++;
			}
			wcout << currentCell->element << endl;
			wcout << endl << "iterations: " << iterations << endl;
			
			return path;
		}
		
		// Find the end cell using a bfs. Print the path taken.
		vector<int> breadthFirstSearch() {
			TreeCell paths = getPaths();
			vector<TreeCell> currentBreadth;
			vector<TreeCell> nextBreadth;
			vector<vector<int>> pathIndeces;
			
			nextBreadth.push_back(paths);
			pathIndeces.push_back({});
			bool hasChildren;
			int bSize;
			
			bool foundEnd = false;
			int endIndex = -1;
			
			int iterations = 0;
			
			while (!foundEnd) {
				// move nextBreadth to currentBreadth
				currentBreadth = {};
				for (int i = 0; i < nextBreadth.size(); i++) {
					currentBreadth.push_back(nextBreadth[i]);
				}
				nextBreadth = {};
				bSize = 0;
				
				/*
				// print the breadth
				for (int i = 0; i < currentBreadth.size(); i++) {
					if (i != 0)
						wcout << ", ";
					wcout << currentBreadth[i].element;
				}
				wcout << endl;
				*/
				
				// for each cell in the current breadth...
				for (int i = 0; i < currentBreadth.size(); i++) {
					
					// If the current cell has the target element, break out of the loop.
					if (currentBreadth[i].element == end) {
						foundEnd = true;
						endIndex = i;
						break;
					}
					
					hasChildren = false;
					// For each of the current cell's children...
					for (int j = 0; j < currentBreadth[i].children.size(); j++) {
						hasChildren = true;
						
						if (j != 0) {
							pathIndeces.insert(pathIndeces.begin() + bSize, pathIndeces[bSize - 1]);
							pathIndeces[bSize].pop_back();
						}
						pathIndeces[bSize].push_back(j);
						nextBreadth.push_back(currentBreadth[i].children[j]);
						bSize = nextBreadth.size();
					}
					
					if (!hasChildren) {
						pathIndeces.erase(pathIndeces.begin() + bSize);
					}
				}
				iterations++;
			}
			wcout << "iterations: " << iterations << endl;
			
			if (endIndex == -1) {
				return {};
			}
			return pathIndeces[endIndex];
		}
		
		// TO-DO: Find the exit using an optimal pathfinding algorithm?
		// - Assume path given minimal obstacles, readjust path when obstacles are encountered.
		
		// Print the maze using box-drawing characters.
		// - NOTE: Each character is an intersection of walls.
		void print() {
			int boxValue;
			bool right;
			
			for (int i = 0; i <= length; i++) {
				for (int j = 0; j <= width * 2; j++) {
					right = false;
					boxValue = 0;
					
					if (i != 0 && verticalWalls[j / 2][i - 1]
					&& ((j != 0 && mask[i - 1][(j-1) / 2]) || (j < width * 2 && mask[i - 1][j / 2]))) {
						
						boxValue += 1;
					}
					
					if (j != 0 && horizontalWalls[i][(j / 2) - 1]
					&& ((i != 0 && mask[i - 1][(j-1) / 2]) || (i < length && mask[i][(j-1) / 2]))) {
						boxValue += 2;
					}
					
					if (j < width * 2 && horizontalWalls[i][j / 2]
					&& ((i != 0 && mask[i - 1][j / 2]) || (i < length && mask[i][j / 2]))) {
						boxValue += 4;
						right = true;
					}
					
					if (i != length && verticalWalls[j / 2][i]
					&& ((j != 0 && mask[i][(j-1) / 2]) || (j < width * 2 && mask[i][j / 2]))) {
						boxValue += 8;
					}
					
					if (j % 2 == 0) {
						switch(boxValue) {
							case 0:
								wcout << L" ";
								break;
							case 1:
								wcout << L"╵";
								break;
							case 2:
								wcout << L"╴";
								break;
							case 3:
								wcout << L"┘";
								break;
							case 4:
								wcout << L"╶";
								break;
							case 5:
								wcout << L"└";
								break;
							case 6:
								wcout << L"─";
								break;
							case 7:
								wcout << L"┴";
								break;
							case 8:
								wcout << L"╷";
								break;
							case 9:
								wcout << L"│";
								break;
							case 10:
								wcout << L"┐";
								break;
							case 11:
								wcout << L"┤";
								break;
							case 12:
								wcout << L"┌";
								break;
							case 13:
								wcout << L"├";
								break;
							case 14:
								wcout << L"┬";
								break;
							case 15:
								wcout << L"┼";
								break;
							default:
								wcout << L" ";
						}
					} else {
						if (right)
							wcout << L"─";
						else
							wcout << " ";
					}
				}
				wcout << endl;
			}
			wcout << endl;
		}
		
		// Print the maze to appear twice as large.
		// - TO-DO: make this function consider the mask
		void printBig() {
			int boxValue;
			bool right;
			bool down;
			for (int i = 0; i <= length * 2; i++) {
				for (int j = 0; j <= width * 4; j++) {
					
					boxValue = 0;
					right = false;
					down = false; 
					
					if (i != 0 && verticalWalls[j / 4][(i / 2) - 1]) {
						boxValue += 1;
					}
					
					if (j != 0 && horizontalWalls[i / 2][(j / 4) - 1]) {
						boxValue += 2;
					}
					
					if (j < length * 4 && horizontalWalls[i / 2][j / 4]) {
						boxValue += 4;
						right = true;
					}
					
					if (i < width * 2 && verticalWalls[j / 4][i / 2]) { 
						boxValue += 8;
						down = true;
					}
					
					if (i % 2 == 0 && j % 4 == 0) {
						switch(boxValue) {
							case 0:
								wcout << L" ";
								break;
							case 1:
								wcout << L"╵";
								break;
							case 2:
								wcout << L"╴";
								break;
							case 3:
								wcout << L"┘";
								break;
							case 4:
								wcout << L"╶";
								break;
							case 5:
								wcout << L"└";
								break;
							case 6:
								wcout << L"─";
								break;
							case 7:
								wcout << L"┴";
								break;
							case 8:
								wcout << L"╷";
								break;
							case 9:
								wcout << L"│";
								break;
							case 10:
								wcout << L"┐";
								break;
							case 11:
								wcout << L"┤";
								break;
							case 12:
								wcout << L"┌";
								break;
							case 13:
								wcout << L"├";
								break;
							case 14:
								wcout << L"┬";
								break;
							case 15:
								wcout << L"┼";
								break;
							default:
								wcout << L" ";
						}	
					}
					
					if (i % 2 == 0 && j % 4 != 0) {
						if (right) {
							wcout << L"─";
						} else {
							wcout << " ";
						}
					}
					
					if (i % 2 != 0 && j % 4 == 0) {
						if (down) {
							wcout << L"│";
						} else {
							wcout << " ";
						}
					}
					
					if (i % 2 != 0 && j % 4 != 0) {
						wcout << " ";
					}
				}
				wcout << endl;
			}
			wcout << endl;
		}
		
		// TO-DO: Print this maze with a solution provided by either search
};

int main(int argc, char *argv[]) {
    _setmode(_fileno(stdout), _O_U16TEXT);
	int seed = time(0);
	// int seed = 1754278747;
	srand(seed);
	wcout << "Seed used: " << seed << endl;
	
	string mask_path = "masks/";
	string ext = ".txt";
	string mask_name = "fish";
	if (argc >= 2) mask_name = argv[1];
	string full_mask_path = mask_path + mask_name + ext;
	
	// Maze maze(20, 50);
	
	vector<vector<bool>> mask = readMaskFromFile(full_mask_path);
	
	Maze maze(mask);
	
	maze.print();
	
	// maze.printBig();
	
	// TreeCell paths = maze.getPaths();
	// paths.print();
	// wcout << endl << endl;
	
	/*
	vector<int> solution1 = maze.depthFirstSearch();
	printVector(solution1);
	wcout << endl << endl;
	
	vector<int> solution2 = maze.breadthFirstSearch();
	printVector(solution2);
	wcout << endl << endl;
	*/
}