using System;


namespace Meta.Numerics.Matrices {

#if FUTURE

    public abstract class MatrixBase<T> {

        protected MatrixBase (int rowCount, int columnCount) {
            this.rowCount = rowCount;
            this.columnCount = columnCount;
        }

        int rowCount, columnCount;

        public virtual int RowCount {
            get {
                return(rowCount);
            }
        }

        public virtual int ColumnCount {
            get {
                return (columnCount);
            }
        }

        protected abstract T GetEntry (int r, int c);

        protected abstract void SetEntry (int r, int c, T value);

        public virtual T this[int r, int c] {
            get {
                return (GetEntry(r, c));
            }
            set {
                SetEntry(r, c, value);
            }
        }

    }

#endif

}
