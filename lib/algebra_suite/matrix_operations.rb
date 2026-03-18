# frozen_string_literal: true

module AlgebraSuite
  # Класс для работы с векторами
  class Vector
    attr_reader :elements

    # Создание вектора
    def initialize(elements)
      raise ArgumentError, "Vector must be an array" unless elements.is_a?(Array)
      raise ArgumentError, "Vector cannot be empty" if elements.empty?
      raise ArgumentError, "Vector elements must be numbers" unless elements.all? { |x| x.is_a?(Numeric) }

      @elements = elements.map(&:to_f)
    end

    # Размер вектора
    def size
      @elements.size
    end

    # Сложение векторов
    def +(other)
      raise ArgumentError, "Vectors must have same size" unless size == other.size

      Vector.new(@elements.each_index.map { |i| @elements[i] + other.elements[i] })
    end

    # Вычитание векторов
    def -(other)
      raise ArgumentError, "Vectors must have same size" unless size == other.size

      Vector.new(@elements.each_index.map { |i| @elements[i] - other.elements[i] })
    end

    # Умножение вектора на число
    def *(scalar)
      raise ArgumentError, "Scalar must be a number" unless scalar.is_a?(Numeric)

      Vector.new(@elements.map { |x| x * scalar })
    end

    # Скалярное произведение
    def dot(other)
      raise ArgumentError, "Vectors must have same size" unless size == other.size

      @elements.each_index.sum { |i| @elements[i] * other.elements[i] }
    end

    # Вернуть обычный массив
    def to_a
      @elements.dup
    end
  end

  # Класс для работы с матрицами
  class Matrix
    attr_reader :rows

    EPSILON = 1e-10

    # Создание матрицы
    def initialize(rows)
      raise ArgumentError, "Matrix must be an array of arrays" unless rows.is_a?(Array) && rows.all? { |r| r.is_a?(Array) }
      raise ArgumentError, "Matrix cannot be empty" if rows.empty? || rows.any?(&:empty?)

      size = rows.first.size
      raise ArgumentError, "All rows must have same length" unless rows.all? { |r| r.size == size }
      raise ArgumentError, "Matrix elements must be numbers" unless rows.flatten.all? { |x| x.is_a?(Numeric) }

      @rows = rows.map { |r| r.map(&:to_f) }
    end

    # Количество строк
    def row_count
      @rows.size
    end

    # Количество столбцов
    def column_count
      @rows.first.size
    end

    # Проверка, что матрица квадратная
    def square?
      row_count == column_count
    end

    # Сложение матриц
    def +(other)
      raise ArgumentError, "Matrices must have same dimensions" unless row_count == other.row_count && column_count == other.column_count

      result = @rows.each_index.map do |i|
        @rows[i].each_index.map do |j|
          @rows[i][j] + other.rows[i][j]
        end
      end

      Matrix.new(result)
    end

    # Вычитание матриц
    def -(other)
      raise ArgumentError, "Matrices must have same dimensions" unless row_count == other.row_count && column_count == other.column_count

      result = @rows.each_index.map do |i|
        @rows[i].each_index.map do |j|
          @rows[i][j] - other.rows[i][j]
        end
      end

      Matrix.new(result)
    end

    # Универсальное умножение
    def *(other)
      if other.is_a?(Numeric)
        Matrix.new(@rows.map { |row| row.map { |x| x * other } })
      elsif other.is_a?(Vector)
        multiply_vector(other)
      elsif other.is_a?(Matrix)
        multiply_matrix(other)
      else
        raise ArgumentError, "Unsupported multiplication"
      end
    end

    # Умножение матрицы на вектор
    def multiply_vector(vector)
      raise ArgumentError, "Matrix columns must match vector size" unless column_count == vector.size

      result = @rows.map do |row|
        row.each_index.sum { |i| row[i] * vector.elements[i] }
      end

      Vector.new(result)
    end

    # Умножение матрицы на матрицу
    def multiply_matrix(other)
      raise ArgumentError, "Matrix sizes are not compatible" unless column_count == other.row_count

      result = Array.new(row_count) { Array.new(other.column_count, 0.0) }

      row_count.times do |i|
        other.column_count.times do |j|
          column_count.times do |k|
            result[i][j] += @rows[i][k] * other.rows[k][j]
          end
        end
      end

      Matrix.new(result)
    end

    # Определитель
    def determinant
      raise ArgumentError, "Determinant is defined only for square matrices" unless square?

      return @rows[0][0] if row_count == 1

      if row_count == 2
        a, b = @rows[0]
        c, d = @rows[1]
        return a * d - b * c
      end

      det = 0.0
      @rows[0].each_index do |j|
        det += ((-1) ** j) * @rows[0][j] * minor(0, j).determinant
      end
      det
    end

    # Минор матрицы
    def minor(i, j)
      Matrix.new(
        @rows.each_with_index.map { |row, r|
          row.each_with_index.map { |x, c| x unless c == j } unless r == i
        }.compact.map(&:compact)
      )
    end

    # Приведение к ступенчатому виду методом Гаусса
    def gaussian_elimination
      matrix = deep_copy(@rows)
      current_row = 0

      column_count.times do |col|
        break if current_row >= row_count

        # Ищем ведущий элемент
        pivot_row = (current_row...row_count).find { |r| matrix[r][col].abs > EPSILON }
        next if pivot_row.nil?

        # Меняем строки местами
        matrix[current_row], matrix[pivot_row] = matrix[pivot_row], matrix[current_row] if pivot_row != current_row

        pivot = matrix[current_row][col]

        # Делим строку на ведущий элемент
        (col...column_count).each do |j|
          matrix[current_row][j] /= pivot
        end

        # Обнуляем элементы ниже ведущего
        ((current_row + 1)...row_count).each do |r|
          factor = matrix[r][col]
          next if factor.abs <= EPSILON

          (col...column_count).each do |j|
            matrix[r][j] -= factor * matrix[current_row][j]
          end
        end

        current_row += 1
      end

      Matrix.new(matrix)
    end

    # Ранг матрицы
    def rank
      matrix = gaussian_elimination.rows

      matrix.count do |row|
        row.any? { |value| value.abs > EPSILON }
      end
    end

    # Обратная матрица через метод Гаусса-Жордана
    def inverse
      raise ArgumentError, "Inverse exists only for square matrices" unless square?

      n = row_count
      left = deep_copy(@rows)
      right = identity_matrix(n)

      n.times do |col|
        # Ищем ведущий элемент
        pivot_row = (col...n).find { |r| left[r][col].abs > EPSILON }
        raise ArgumentError, "Matrix is singular and has no inverse" if pivot_row.nil?

        # Меняем строки местами
        if pivot_row != col
          left[col], left[pivot_row] = left[pivot_row], left[col]
          right[col], right[pivot_row] = right[pivot_row], right[col]
        end

        pivot = left[col][col]

        # Нормализуем строку
        n.times do |j|
          left[col][j] /= pivot
          right[col][j] /= pivot
        end

        # Обнуляем остальные элементы в столбце
        n.times do |r|
          next if r == col

          factor = left[r][col]
          next if factor.abs <= EPSILON

          n.times do |j|
            left[r][j] -= factor * left[col][j]
            right[r][j] -= factor * right[col][j]
          end
        end
      end

      Matrix.new(right)
    end

    # Решение системы Ax = b методом Гаусса
    def solve(vector)
      raise ArgumentError, "Matrix must be square" unless square?
      raise ArgumentError, "Expected a Vector" unless vector.is_a?(Vector)
      raise ArgumentError, "Vector size must match matrix size" unless vector.size == row_count

      n = row_count
      augmented = @rows.each_index.map { |i| @rows[i] + [vector.elements[i]] }

      n.times do |col|
        # Ищем ведущий элемент
        pivot_row = (col...n).find { |r| augmented[r][col].abs > EPSILON }
        raise ArgumentError, "System has no unique solution" if pivot_row.nil?

        # Меняем строки местами
        augmented[col], augmented[pivot_row] = augmented[pivot_row], augmented[col] if pivot_row != col

        pivot = augmented[col][col]

        # Нормализуем строку
        (col..n).each do |j|
          augmented[col][j] /= pivot
        end

        # Обнуляем остальные элементы столбца
        n.times do |r|
          next if r == col

          factor = augmented[r][col]
          next if factor.abs <= EPSILON

          (col..n).each do |j|
            augmented[r][j] -= factor * augmented[col][j]
          end
        end
      end

      Vector.new(n.times.map { |i| augmented[i][n] })
    end

    # Вернуть обычный массив
    def to_a
      @rows.map(&:dup)
    end

    private

    def deep_copy(array)
      array.map(&:dup)
    end

    def identity_matrix(size)
      Array.new(size) do |i|
        Array.new(size) { |j| i == j ? 1.0 : 0.0 }
      end
    end
  end
end