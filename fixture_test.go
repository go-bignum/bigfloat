// Copyright 2026 the go-bignum/bigfloat authors. All rights reserved.
// Use of this source code is governed by a BSD-style licence that can be
// found in the LICENSE file.

package bigfloat_test

import (
	"math/big"
	"path/filepath"
	"slices"
	"strings"
	"testing"

	"github.com/go-bignum/bigfloat/oracle"
)

// TestFixtures compares every function, in every rounding mode, with the
// golden results that cmd/oracle generated from mpmath: value, sign of
// zero and accuracy must all be equal — correct rounding admits exactly one
// answer.
func TestFixtures(t *testing.T) {
	files, err := filepath.Glob("testdata/*.json.gz")
	if err != nil || len(files) == 0 {
		t.Fatalf("no fixtures: %v", err)
	}
	names := map[string]bool{}
	for _, path := range files {
		f, err := oracle.LoadFixture(path)
		if err != nil {
			t.Fatal(err)
		}
		names[f.Fn] = true
		t.Run(f.Fn, func(t *testing.T) {
			t.Parallel()
			bad := 0
			for i := range f.Cases {
				c := &f.Cases[i]
				if testing.Short() && c.Prec > 200 {
					continue
				}
				for _, msg := range c.Check(f.Fn, func(z *big.Float, a ...*big.Float) *big.Float { return oracle.Call(f.Fn, z, a...) }) {
					if bad++; bad <= 10 {
						t.Error(msg)
					}
				}
			}
			if bad > 10 {
				t.Errorf("... %d mismatches in all", bad)
			}
		})
	}
	var missing []string
	for fn := range oracle.Functions {
		if !names[fn] {
			missing = append(missing, fn)
		}
	}
	slices.Sort(missing)
	if len(missing) > 0 {
		t.Errorf("functions without fixtures: %s", strings.Join(missing, ", "))
	}
}
