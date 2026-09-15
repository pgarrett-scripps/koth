//! Background-thread prefetching for streaming spectrum sources.
//!
//! mzML reading is decompression + XML parsing bound (~half of a typical
//! run is spent inside `miniz_oxide::inflate` and `quick_xml`), and it runs
//! on the same thread that consumes the spectra for hill detection. Those two
//! costs are independent, so serializing them wastes wall-clock equal to the
//! smaller of the two.
//!
//! [`prefetch`] moves the source iterator onto a dedicated reader thread that
//! decodes spectra into a bounded channel while the caller runs `process_scan`
//! on the previously-decoded one. Spectra are yielded in the exact same order
//! (single producer → single consumer FIFO), so hill detection stays
//! byte-for-byte deterministic — this only overlaps the two stages in time.

use std::sync::mpsc::{sync_channel, Receiver};
use std::thread::JoinHandle;

/// Run a fallible producer on a background thread. Returning `false` from the
/// supplied emitter means the consumer has stopped; producers must stop too.
/// Initialization and decode errors are delivered through the same bounded
/// channel as the data, so they cannot silently truncate a successful result.
pub fn try_prefetch<T, E, F>(producer: F, capacity: usize) -> Prefetch<Result<T, E>>
where
    T: Send + 'static,
    E: Send + 'static,
    F: FnOnce(&mut dyn FnMut(T) -> bool) -> Result<(), E> + Send + 'static,
{
    let (tx, rx) = sync_channel(capacity);
    let handle = std::thread::Builder::new()
        .name("koth-reader".into())
        .spawn(move || {
            if let Err(error) = producer(&mut |item| tx.send(Ok(item)).is_ok()) {
                let _ = tx.send(Err(error));
            }
        })
        .expect("failed to spawn reader thread");
    Prefetch {
        rx: Some(rx),
        handle: Some(handle),
    }
}

/// Run `source` on a dedicated background thread, buffering up to `capacity`
/// items in a bounded channel. The returned iterator yields the same items in
/// the same order.
///
/// The channel is bounded so a fast reader can't run away with memory: once
/// `capacity` spectra are queued the reader blocks until the consumer drains
/// one. A producer panic is re-raised on the consumer thread (see
/// [`Prefetch::next`]) so parse failures still fail loud.
pub fn prefetch<I>(source: I, capacity: usize) -> Prefetch<I::Item>
where
    I: Iterator + Send + 'static,
    I::Item: Send + 'static,
{
    let (tx, rx) = sync_channel(capacity);
    let handle = std::thread::Builder::new()
        .name("koth-prefetch".into())
        .spawn(move || {
            for item in source {
                // Err means the consumer was dropped — stop reading early.
                if tx.send(item).is_err() {
                    break;
                }
            }
        })
        .expect("failed to spawn prefetch thread");
    Prefetch {
        rx: Some(rx),
        handle: Some(handle),
    }
}

/// Iterator over spectra produced by a background reader thread.
pub struct Prefetch<T> {
    rx: Option<Receiver<T>>,
    handle: Option<JoinHandle<()>>,
}

impl<T> Prefetch<T> {
    /// Drop the receiver (signalling the producer to stop) and join the reader
    /// thread, re-raising its panic if it panicked.
    fn finish(&mut self) {
        // Drop the receiver first so a producer blocked on a full channel
        // unblocks (its `send` returns Err); otherwise the join would deadlock.
        self.rx = None;
        if let Some(h) = self.handle.take() {
            if let Err(panic) = h.join() {
                std::panic::resume_unwind(panic);
            }
        }
    }
}

impl<T> Iterator for Prefetch<T> {
    type Item = T;

    fn next(&mut self) -> Option<T> {
        match self.rx.as_ref().and_then(|rx| rx.recv().ok()) {
            Some(item) => Some(item),
            None => {
                // Sender dropped → source exhausted (or the reader panicked).
                // Join so a producer panic surfaces here instead of silently
                // truncating the spectrum stream.
                self.finish();
                None
            }
        }
    }
}

impl<T> Drop for Prefetch<T> {
    fn drop(&mut self) {
        // Best-effort cleanup if the consumer abandons iteration early. Ignore
        // a producer panic here — re-raising during an unwind would abort.
        self.rx = None;
        if let Some(h) = self.handle.take() {
            let _ = h.join();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::{
        atomic::{AtomicUsize, Ordering},
        Arc,
    };
    use std::time::Duration;

    #[test]
    fn fallible_producer_is_bounded_and_cancels_on_drop() {
        let count = Arc::new(AtomicUsize::new(0));
        let produced = count.clone();
        let (ready_tx, ready_rx) = std::sync::mpsc::channel();
        let stream = try_prefetch::<usize, (), _>(
            move |emit| {
                for i in 0..1000 {
                    produced.fetch_add(1, Ordering::SeqCst);
                    if i == 4 {
                        ready_tx.send(()).unwrap();
                    }
                    if !emit(i) {
                        break;
                    }
                }
                Ok(())
            },
            4,
        );
        ready_rx.recv_timeout(Duration::from_secs(5)).unwrap();
        // Four queued values and one blocked send, independent of run length.
        assert_eq!(count.load(Ordering::SeqCst), 5);
        drop(stream); // must unblock the send and join the producer
        assert_eq!(count.load(Ordering::SeqCst), 5);
    }

    #[test]
    fn fallible_producer_delivers_ordered_values_then_error() {
        let mut stream = try_prefetch(
            |emit| {
                assert!(emit(1));
                assert!(emit(2));
                Err("decode failure")
            },
            1,
        );
        assert_eq!(stream.next(), Some(Ok(1)));
        assert_eq!(stream.next(), Some(Ok(2)));
        assert_eq!(stream.next(), Some(Err("decode failure")));
        assert_eq!(stream.next(), None);
        assert_eq!(stream.next(), None);
    }

    #[test]
    fn producer_panic_reaches_consumer() {
        let mut stream = try_prefetch::<(), (), _>(|_| panic!("reader panic"), 1);
        assert!(std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| stream.next())).is_err());
    }
}
